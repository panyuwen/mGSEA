import numpy as np
import pandas
import argparse
from concurrent.futures import ProcessPoolExecutor

## geneset sampling, to generate the NULL distribution
def genesetperm(genelist,perm,setsize):
    ## sampling
    sample = pandas.DataFrame()
    for p in range(perm):
        sample['p'+str(p)] = list(np.random.choice(genelist,size=setsize,replace=False))

    return sample

## for permutation
def scoresummarize(rnk,qsize,setgene,marker):
    qbed = quantilescore(rnk,qsize,setgene)
    qbed.index = qbed['qstart'].apply(str) +'-'+ qbed['qend'].apply(str)

    ## summarize
    summary = pandas.DataFrame()
    summary.loc[0,'ID'] = marker

    for s in ['hitscore','misscore']:
        for index in [1,2,3]:
            scoretype = s + '_cum'
            q = 'quarter' + str(index)
            quantile = index * 0.25

            summary.loc[0,s+'_'+q+'_region'] = list(qbed[qbed[scoretype]>=quantile].index)[0]
    
    summary.loc[0,'netscore_max'] = qbed['netscore_cum'].max()
    summary.loc[0,'netscore_min'] = qbed['netscore_cum'].min()

    return summary

## calculate cumulative scores for each quantile region
## for real data
def quantilescore(qrnk,qsize,setgene):
    rnk = qrnk.copy()
    rnk['geneset'] = np.isin(rnk['gene'].values, setgene) - 0

    ## quantile region
    qbed = pandas.DataFrame(columns=['qstart','qend'])
    qbed['qstart'] = rnk['quantile'].unique()
    qbed['qend'] = qbed['qstart'].values + qsize
    qbed.index = qbed['qstart'].values

    ## total genes and set genes in each quantile region
    qtot = pandas.value_counts(rnk['quantile'].values)
    qset = rnk[['quantile','geneset']].groupby('quantile').sum()

    qbed = pandas.concat([qbed, qtot, qset],axis=1)
    qbed.columns = ['qstart','qend','gene_tot','gene_set']
    qbed.sort_values(by=['qstart'],ascending=False,inplace=True)

    qbed['gene_mis'] = qbed['gene_tot'] - qbed['gene_set']

    ## weight score, cumulative score
    rnksize = rnk.shape[0]; setsize = len(setgene)
    hitscore = 1.0 / setsize; misscore = 1.0 / (rnksize - setsize)

    #qbed['hitscore'] = qbed['gene_set'] * hitscore; qbed['misscore'] = qbed['gene_mis'] * misscore
    #qbed['netscore'] = qbed['hitscore'] - qbed['misscore']
    qbed['hitscore_cum'] = (qbed['gene_set'] * hitscore).cumsum()
    qbed['misscore_cum'] = (qbed['gene_mis'] * misscore).cumsum()
    qbed['netscore_cum'] = qbed['hitscore_cum'] - qbed['misscore_cum']

    return qbed

def qgsea(rnk,qsize,perm,gsetsample,nthread,output,setname,setdesp,*setgene):
    genelist = list(rnk['gene']); setgene = list(set(genelist) & set(setgene))
    rnksize = len(genelist); setsize = len(setgene)

    ## cumulative score for observed data
    qscore = quantilescore(rnk,qsize,setgene)
    qscore.to_csv(output+'.'+setname+'.cumulative_score.txt.gz',sep='\t',index=None,compression='gzip')

    ## summary of permuation
    with ProcessPoolExecutor(max_workers = nthread) as pool:
        qsummary = pandas.concat(list(pool.map(scoresummarize, [rnk]*(perm+1), [qsize]*(perm+1), [setgene]+[list(gsetsample['p'+str(p)]) for p in range(perm)],['obs']+['p'+str(p+1) for p in range(perm)])),ignore_index=True)
    
    qsummary.index = qsummary['ID']
    qsummary['netscore_sign'] = qsummary.apply(lambda x: 'pos' if (abs(x['netscore_max'])>abs(x['netscore_min'])) else 'neg',axis=1)
    
    qsummary.to_csv(output+'.'+setname+'.enrichment_score_quarter.txt.gz',sep='\t',index=None,compression='gzip')
    
    return qscore, qsummary

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qrank", type=str, required=True, \
                        help="<gene> <quantile value>, tab seperated, no header, no duplicated genes. better to keep only genes in the background gene set, e.g., KEGG genes.")
    parser.add_argument("--gset", type=str, required=True, \
                        help="Priori target gene list, all of the genes should be in the ranked genelist, gmt format: <gene set name> <description> - <gene list> -")
    parser.add_argument("--qsize", type=float, required=True, \
                        help="Window size of the non-overlapping quantile regions, from 1.0 to 0.0, e.g., 0.005")
    parser.add_argument("--perm", type=int, required=False, default=2000, \
                        help="Times of permutation, Default: 2000")
    parser.add_argument("--threads", type=int, required=False, default=50, \
                        help="Number of threads, Default: 50")
    parser.add_argument("--out", type=str, required=True, \
                        help="prefix name for output file.")
    args = parser.parse_args()

    ## genelist ranked according to their quantiles
    qrank = pandas.read_csv(args.qrank,sep='\t',header=None,names=['gene','quantile'])
    qrank.dropna(inplace=True)
    qrank.sort_values(by=['quantile'],ascending=False,inplace=True)
    qrank.drop_duplicates(['gene'],inplace=True)
    allgenes = list(qrank['gene'])
    ## group into quantile regions
    qrank['quantile'] = (qrank['quantile'] / args.qsize).astype(int) * args.qsize

    ## group genesets according to their gene sizes
    with open(args.gset) as gset:
        genesetsize = {}
        line = gset.readline()
        while line:
            tmpsize = len(set(line.strip().split('\t')[2:]) & set(allgenes))
            if tmpsize >= 5:
                if tmpsize not in list(genesetsize.keys()):
                    genesetsize[tmpsize] = [line]
                else:
                    genesetsize[tmpsize].append(line)
            else:
                pass
            line = gset.readline()

    if genesetsize == {}:
        print('no enough set genes.')
        exit()
    else:
        pass

    statistic = pandas.DataFrame(); n = 1
    #with open(args.out+'.enrichment_Pvalue.txt','w') as f:
    #    f.write("geneset\tdescription\thitscore_empirical.Pvalue\tmisscore_empirical.Pvalue\tnetscore_1tail.Pvalue\tnetscore_2tail.Pvalue\tleading.edge_gene\n")

    for setsize in list(genesetsize.keys()):
        genesetsample = genesetperm(list(qrank['gene']),args.perm,setsize)
        ## geneset of the same size
        for line in genesetsize[setsize]:
            setname,setdesp = line.strip().split('\t')[:2]
            setgenelist = line.strip().split('\t')[2:]
            ## summary statistic for observed data and permutated data
            ## the quantile region that firstly achieve the cumulative hit/miss score of 0.25, 0.5, and 0.75
            ## enrichment score
            qscore, qsummary = qgsea(qrank,args.qsize,args.perm,genesetsample,args.threads,args.out,*(line.strip().split('\t')))
            ## estimate the significance based on the summary statistic
            statistic.loc[n,'geneset'] = setname; statistic.loc[n,'description'] = setdesp
            statistic.loc[n,'geneset_size'] = setsize
            ## cumulative hit/miss score
            #for s in ['hitscore','misscore']:
            #    qsummary[s+'_quarter2_start'] = qsummary[s+'_quarter2_region'].apply(lambda x: float(x.split('-')[0]))
            #    statistic.loc[n,s+'_empirical.Pvalue'] = qsummary[qsummary[s+'_quarter2_start']>=qsummary.loc['obs',s+'_quarter2_start']].shape[0]*1.0 / (qsummary.shape[0] + 1.0)
            ## enrichment score
            ## one-tail P-value suggested
            obs_net_max = qsummary.loc['obs','netscore_max']
            obs_net_min = qsummary.loc['obs','netscore_min']

            pos_sign = qsummary[qsummary['netscore_sign']=='pos'].copy()
            pos_sign_signif = pos_sign[pos_sign['netscore_max']>=obs_net_max].copy()
            neg_sign = qsummary[qsummary['netscore_sign']=='neg'].copy()
            neg_sign_signif = neg_sign[neg_sign['netscore_min']<=obs_net_min].copy()

            if qsummary.loc['obs','netscore_sign'] == 'pos':
                qscore.index = qscore['qstart'].values
                statistic.loc[n,'netscore_1tail.Pvalue'] = pos_sign_signif.shape[0]*1.0 / (pos_sign.shape[0] + 1.0)
                peak = qscore['netscore_cum'].idxmax()
                leadgene = ';'.join(list(qrank[(qrank['quantile']>=peak) & (qrank['gene'].isin(setgenelist))]['gene']))
            else:
                qscore.index = qscore['qend'].values
                statistic.loc[n,'netscore_1tail.Pvalue'] = neg_sign_signif.shape[0]*1.0 / (neg_sign.shape[0] + 1.0)
                peak = qscore['netscore_cum'].idxmin()
                leadgene = ';'.join(list(qrank[(qrank['quantile']<peak) & (qrank['gene'].isin(setgenelist))]['gene']))

            statistic.loc[n,'netscore_2tail.Pvalue'] = min(pos_sign_signif.shape[0], neg_sign_signif.shape[0]) * 2.0 / (qsummary.shape[0] + 1.0)
            statistic.loc[n,'leading.edge_gene'] = leadgene

            #with open(args.out+'.enrichment_Pvalue.txt','a') as f:
            #    f.write('\t'.join(statistic.loc[n].apply(str))+'\n')

            n += 1

    statistic.sort_values(by=['netscore_1tail.Pvalue'],ascending=True,inplace=True)
    statistic.to_csv(args.out+'.enrichment_Pvalue.txt',sep='\t',index=None)

if __name__ == '__main__':
    main()
