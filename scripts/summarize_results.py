import pandas as pd
import numpy as np
from collections import Counter
from Bio import Phylo

def process_blast_recomb(name):
    """Extract results from sliding window BLAST

    This function parses the results from the slidding window BLAST
    ('recblast_{id_of_the_fasta_sequence}.txt'), filters the top results,
    divide the results by bins (of 50 nucleotides), and outputs the majority
    rule result for r«each been. Finaly an output for the given sequence is
    returned, if there is multiple results they are separated by '/'.


	Args:
        name (txt): Sliding window BLAST result.

	Returns:
        List with two items: name and processed result.
	"""
    name_out = name[:-4].replace('blast/recblast_', '')

    try:
        df = pd.read_csv(name, header=None)
        df = df[df[2] < 0.1e-105].copy()
    except:
        return [name_out, [np.NaN]]
    
    dif_splits = df[0].unique()
    as_array = []

    for split in dif_splits:
        position = int(split.split('_')[-1])

        result = list(int(position / 50) * '-')
        
        split_df = df[df[0] == split]

        split_df = split_df.sort_values(by=[3], ascending=False)

        try:
            best_score = split_df[3].values[0]
            refs_best_score = split_df[split_df[3] == best_score][1].values
            subs_best = list(set([x.split('-')[0] for x in refs_best_score]))

            if len(subs_best) == 1:
                 subtype = subs_best[0]
            else:
                subtype = '-'

            result += [subtype] * 8
            
        except:
            pass

        bins_number = (int(dif_splits[-1].split('_')[-1]) + 400) / 50
        result += list('-' * int(bins_number - len(result)))

        as_array.append(np.array(result)) 

    most_comon_in_array = [Counter([i for i in x if i != '-']).most_common(1) for x in np.array(as_array).T]
    
    pass_res = [x[0] for x in most_comon_in_array if len(x) == 1]
    final = list(set([x[0] for x in pass_res if x[1] > 4]))
    
    return [name_out, ["/".join(final)]]



def process_trees(tree):
    """Extract results from the phylogenetic inference

    This function parses the results from the phylogenetic trees created, roots
    the trees on the outgrouop ('CONSENSUS_CPZ'), and evaluated if the target
    sequence is or not in a monophyletic clade with reference of only one
    subtype/circulating recombinat form (crf). If the previously stated conditions
    occur a subtype/crf is outputed, otherwise np.nan is outputed.

	Args:
        tree (nwk): Phylogenetic tree in nwk format.

	Returns:
        List with two items: name and processed result.
	"""
    name_target = tree[:-4].replace('trees/all_', '').replace('trees/pure_', '').replace('trees/recomb_', '')

    t = Phylo.read(tree, 'newick')
    t.root_with_outgroup('CONSENSUS_CPZ')
    t.ladderize()

    nodes_with_target = t.get_path(name_target)

    result = []
    for node in nodes_with_target:
        in_node = list(set([i.name.split('-')[0] for i in node.get_terminals() if i.name != name_target]))
        if len(in_node) == 1:
            result = [name_target, in_node[0], node.confidence]
            break
        else:
            pass
        
    if result == []:
        result = [name_target, np.NaN, 0]
    else:
        pass
    
    return result


def get_clossest_blast(blast_res):
    """Extract closser reference from the BLAST

    This function parses the results from the BLAST, turn it into a pandas
    dataframe, verifies how many diferent subtype sequences have the top BLAST
    score, if more than one the output is np.nan, if only one the output is that
    subtype/CRF. 

	Args:
        blast_res (txt): BLAST result. Is csv.

	Returns:
        List with two items: name and processed result.
	"""
    name_out = blast_res.replace('blast/blast_', '').replace('.txt', '')

    df = pd.read_csv(blast_res)
    df.columns = [0,1,2,3]

    filter_df = df.sort_values(by=[3], ascending=False).copy()

    best_score = filter_df[3].values[0]

    refs_best_score = filter_df[filter_df[3] == best_score][1].values
    subs_best = list(set([x.split('-')[0] for x in refs_best_score]))

    if len(subs_best) == 1:
        return [name_out, subs_best[0]]
    else:
        return [name_out, np.NaN]

def make_decision(idx, df):
    """Create final result based on all analysis 

    This function consists in a series of if statments. If statement can be seen
    as a 'rule' with requirements that need to be meet to the final SNAPPy
    result to be created.  

	Args:
        idx (str): Internal SNAPPy id.
        df (dataframe): Tabular like file with the outputs from all the analysis
        performed.

	Returns:
        List with two items: rule used and final SNAPPy output.
	"""

    to_process = list(df.loc[idx])

    # all methods agree
    ## rule_p1: no recomb, tree all equal tree pure, recomb equal tree all, tree all equal closser
    if ((to_process[8] == 0) & (str(to_process[1]) != 'nan') &
        (to_process[1] == to_process[3]) & (to_process[0] == to_process[1]) &
        (to_process[1] == to_process[7])):
        return ['rule_p1', to_process[1]]
    ## rule_c1: all trees and recomb trees and closser ref agree plus recomb is simple 
    elif ((str(to_process[1]) != 'nan') & (to_process[1] == to_process[5]) &
          (to_process[2] >= 0.7) & (to_process[6] >= 0.7) &
          (to_process[8] == 1) & (to_process[1] == to_process[7]) &
          (str(to_process[1]) != 'nan')):
        return ['rule_c1', to_process[1]]
    
    # both trees plus 1 method agree
    ## rule_p2: tree pure agrees with tree all and recomb
    elif ((str(to_process[3]) != 'nan') & (to_process[3] == to_process[1]) &
        (to_process[4] >=0.7) & (to_process[2] >=0.7) &
        (to_process[3] == to_process[0])):
        return ['rule_p2', to_process[3]]
    ## rule_p3: tree pure agrees with tree all and closser
    elif ((str(to_process[3]) != 'nan') & (to_process[3] == to_process[1]) &
        (to_process[4] >=0.7) & (to_process[2] >=0.7) &
        (to_process[3] == to_process[7])):
        return ['rule_p3', to_process[3]]  
    ## rule_c2: tree recomb agrees with tree all and closser and there is recomb
    elif ((str(to_process[5]) != 'nan') & (to_process[5] == to_process[1]) &
        (to_process[6] >=0.7) & (to_process[2] >=0.7) &
        (to_process[5] == to_process[7])):
        return ['rule_c2', to_process[5]]

    # one tree plus recomb and closser
    ## rule_p4: tree pure agrees with recomb and closser
    elif ((str(to_process[3]) != 'nan') & (to_process[4] >=0.9) & 
        (to_process[3] == to_process[0]) & (to_process[3] == to_process[7])):
        return ['rule_p4', to_process[3]]
    ## rule_c3: tree recomb agrees with closser and recomb is simple 
    elif ((str(to_process[5]) != 'nan') & (to_process[6] >=0.9) & 
        (to_process[8] == 1) & (to_process[5] == to_process[7])):
        return ['rule_c3', to_process[5]]
    ## rule_b1: tree all agrees with recomb and closser
    elif ((str(to_process[1]) != 'nan') & (to_process[2] >=0.9) & 
        (to_process[1] == to_process[0]) & (to_process[1] == to_process[7])):
        return ['rule_b1', to_process[1]]
    
    # ecomb gives complex
    ## rules_c4: tree all agrees tree recomb, and their result is a crf
    elif ((to_process[8] == 2)):
        if ((to_process[1] == to_process[5]) & (to_process[2] >= 0.7) &
            (to_process[6] >= 0.7) & ('_' in str(to_process[1]))):
            return ['rule_c4', to_process[1]]
        ## rules_p5: tree all agrees tree pure, and closser, great support for 1 tree
        elif ((to_process[1] == to_process[3]) & (to_process[1] == to_process[7]) &
              ((to_process[2] >= 0.9) | (to_process[4] >=0.9))):
            return ['rule_p5', to_process[1]]
        ## rules_c5: tree all agrees tree recomb, and closser, and trees give crf  
        elif ((to_process[1] == to_process[5]) & ('_' in str(to_process[1])) &
              (to_process[1] == to_process[7])):
            return ['rule_c5', to_process[1]]
        ## rules_p6: tree all agrees tree pure, and closser
        elif ((to_process[1] == to_process[3]) & (to_process[1] == to_process[7])):
            return ['rule_p6', to_process[1]]
        ## rules_u1: remaining cases are an complex URF       
        else:     
            return ['rule_u1', 'URF_CPX']
    
    # recomb gives simple
    ## rules_c6: tree all agrees tree recomb, and their result is a crf
    elif ((to_process[8] == 1)):
        if ((to_process[1] == to_process[5]) & (to_process[2] >= 0.7) &
            (to_process[6] >= 0.7) & ('_' in str(to_process[1]))):
            return ['rule_c6', to_process[1]]
        ## rules_p7: tree all agrees tree pure, and closser, great support for 1 tree
        elif ((to_process[1] == to_process[3]) & (to_process[1] == to_process[7]) &
              ((to_process[2] >= 0.9) | (to_process[4] >=0.9))):
            return ['rule_p7', to_process[1]]
        ## rules_c7: tree all agrees tree recomb, and closser, and trees give crf     
        elif ((to_process[1] == to_process[5]) & ('_' in str(to_process[1])) &
              (to_process[1] == to_process[7])):
            return ['rule_c7', to_process[1]]
        ## rules_p8: tree all agrees tree pure, and closser    
        elif ((to_process[1] == to_process[3]) & (to_process[1] == to_process[7])):
            return ['rule_p8', to_process[1]]
        ## rules_u1: remaining cases are an URF     
        else:
            return ['rule_u2', f'URF_{"".join([str(x)[:2] for x in sorted(str(to_process[0]).split("/"))])}']
                
    # no evidence of recomb
    ## rule_p9: pure and all trees agree
    elif ((to_process[1] == to_process[3]) &
         (to_process[4] >=0.7) & (to_process[2] >=0.7)):
         return ['rule_p9', to_process[1]]
      
    # final, deal with problems of missing data
    ## rule_f1: if recomb res mssing output closser sed result
    elif ((str(to_process[0]) == 'nan') & (str(to_process[7]) != 'nan')):
        return ['rule_f1', to_process[7]]
    ## rule_f2: if recomb res and closser outputs misisng and trees agree give trees result
    elif ((str(to_process[0]) == 'nan') & (str(to_process[7]) == 'nan')):
        if ((to_process[1] == to_process[3]) & (str(to_process[1]) != 'nan')):
            return ['rule_f2', to_process[1]]
        ## rule_f3: if recomb res and closser outputs misisng and trees agree give trees result
        elif ((to_process[1] == to_process[5]) & (str(to_process[1]) != 'nan')):
            return ['rule_f3', to_process[1]] 
        ## rule_f4: else return impossible to dtermine
        else:
            return ['rule_f4', 'impossible_to_determine'] 
                    
    ## rule_f5: give what is ouputed be recomb test, there is no recomb
    else:
        return ['rule_f5', to_process[0]]


def has_recomb(str_recomb):
    n_recombs = len(str(str_recomb).split('/'))
    if n_recombs == 1:
        return 0
    elif n_recombs == 2:
        return 1
    elif n_recombs > 2:
        return 2
    else:
        print('ERRRRROOR')

if __name__ == '__main__':
    input_list = list(snakemake.input)
    keys = list(snakemake.params.k)
    ids = list(snakemake.params.i)

    all_trees_inputs = [x for x in input_list if (x[-4:]  == '.nwk') & (x[:10] == 'trees/all_')]
    pure_trees_inputs = [x for x in input_list if (x[-4:]  == '.nwk') & (x[:11] == 'trees/pure_')]
    recomb_trees_inputs = [x for x in input_list if (x[-4:]  == '.nwk') & (x[:13] == 'trees/recomb_')]
    blast_c = [x for x in input_list if (x[-4:]  == '.txt') & (x[:12] == 'blast/blast_')]
    blast_inputs = [x for x in input_list if (x[-4:]  == '.txt') & (x[:15] == 'blast/recblast_')]

    results = {}
    for pos, key in enumerate(keys):
        results[key] = [ids[pos]]

    for blast_res in blast_inputs:
        output = process_blast_recomb(blast_res)
        results[output[0]] += output[1]

    for tree in all_trees_inputs:
        output = process_trees(tree)
        results[output[0]] += output[1:]
    
    for tree in pure_trees_inputs:
        output = process_trees(tree)
        results[output[0]] += output[1:]

    for tree in recomb_trees_inputs:
        output = process_trees(tree)
        results[output[0]] += output[1:]

    for blast in blast_c:
        output = get_clossest_blast(blast)
        results[output[0]] += [output[1]]

    df_report = pd.DataFrame.from_dict(results, orient='index')

    df_report.columns = ['id', 'recomb_result', 'node_all_refs', 's_node_all_refs',
                      'node_pure_refs', 's_node_pure_refs', 'node_recomb_refs',
                      's_node_recomb_refs', 'closser_ref']


    
    to_make_decision = df_report.set_index(['id']).copy()
    to_make_decision['has_recomb'] = to_make_decision['recomb_result'].apply(lambda x: has_recomb(x))

    result_dict = {}
    for idx in list(to_make_decision.index):
        result_dict[idx] = make_decision(idx, to_make_decision)

    my_res = pd.DataFrame.from_dict(result_dict, orient='index')
    my_res.columns = ['rule', 'my_result']
    my_res['my_result'] = my_res['my_result'].str.upper()
    my_res['my_result'] = my_res['my_result'].str.replace('32_06A6', '32_06A1')
    my_res['my_result'] = my_res['my_result'].str.replace('URF_A1A2', 'A')
    my_res['my_result'] = my_res['my_result'].str.replace('URF_F1F2', 'F')

    my_res.reset_index(inplace=True)
    df_report.reset_index(inplace=True)

    df_report['rule'] = my_res['rule']
    df_report['result'] = my_res['my_result']

    df_report.columns = ['id', 'name', 'recomb_result', 'node_all_refs',
                         's_node_all_refs', 'node_pure_refs', 's_node_pure_refs',
                         'node_recomb_refs', 's_node_recomb_refs', 'closser_ref',
                         'rule', 'result']

    df_report[['id', 'name', 'result', 'recomb_result', 'node_all_refs',
              's_node_all_refs', 'node_pure_refs', 's_node_pure_refs',
              'node_recomb_refs', 's_node_recomb_refs', 'closser_ref',
              'rule']].to_csv('report_subtype_results.csv', index=None)

    my_res['id'] = df_report['id']
    my_res.columns = ['name', 'rule', 'result', 'id']

    my_res[['id', 'name', 'result']].to_csv('subtype_results.csv', index=None)
