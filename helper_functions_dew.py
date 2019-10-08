
import pickle
import os
import scipy.sparse
import numpy as np
import pandas as pd
import scanpy.api as sc
import matplotlib as mpl
import matplotlib.pyplot as plt


# LOADING DATA

def load_inDrops_V3(library_names, input_path):
    '''
    Imports inDrops V3 data files.  The first time this function is executed, it will load
    counts matrices, gene names, cell names, and cell barcode sequences from original tsv and pickle
    files, respectively.  Fast-loading versions of these objects (e.g. *.npz) will be saved for 
    future calls to this function.
    The returned dictionary object D includes the following entries: 
    'E', meta', 'gene_names', 'cell_names', 'cell_bc_seqs'
    '''

    # Create a dictionary to hold data
    D = {}
    for j, s in enumerate(library_names):
        D[s] = {}

    # Load counts data, metadata, & convert to AnnData objects
    for s in library_names:
        print('_________________', s)

        # First attempt to load matrix data from preprocessed files (fast)
        if os.path.isfile(input_path + s + '/' + s + '.raw_counts.unfiltered.npz'):
            print('Loading from npz file')
            E = scipy.sparse.load_npz(
                input_path + s + '/' + s + '.raw_counts.unfiltered.npz')
            gene_names = np.loadtxt(
                fname=input_path + s + '/gene_names.txt', dtype='str')
            cell_names = np.loadtxt(
                fname=input_path + s + '/cell_names.txt', dtype='str')
            cell_bc_seqs = np.loadtxt(
                fname=input_path + s + '/cell_bc_seqs.txt', dtype='str')

        # Otherwise, load and preprocess from the original text files (slow)
        else:
            print('Loading from text file')
            counts_mat = pd.read_csv(
                input_path + s + '/' + s + '.counts.tsv.gz', sep='\t', index_col=0)
            E = scipy.sparse.coo_matrix(np.asmatrix(counts_mat.values)).tocsc()
            cell_names = counts_mat.index
            gene_names = counts_mat.columns

            # Load the barcode dictionary pickle file, format as keys=bcodes; values=sequences
            f = open(input_path + s + '/abundant_barcodes.pickle', 'rb')
            bc_dict = pickle.load(f)
            f.close()
            bcd_dict = {bc_dict[bc][0]: bc for bc in bc_dict}

            # Get barcode sequences corresponding to each cell index
            bcd_seqs = []
            for cname in counts_mat.index:
                bcd_seqs.append(s + '_' + bcd_dict.get(cname))
            cell_bc_seqs = bcd_seqs

            # Save fast files for next time
            scipy.sparse.save_npz(input_path + s + '/' +
                                  s + '.raw_counts.unfiltered.npz', E)
            np.savetxt(input_path + s + '/gene_names.txt',
                       counts_mat.columns, fmt='%s')
            np.savetxt(input_path + s + '/cell_names.txt',
                       counts_mat.index, fmt='%s')
            np.savetxt(input_path + s + '/cell_bc_seqs.txt',
                       bcd_seqs, fmt='%s')

        # Print matrix dimensions to screen
        print(E.shape, '\n')

        # Convert to ScanPy AnnData objects
        D[s]['adata'] = sc.AnnData(E)
        D[s]['adata'].obs['n_counts'] = D[s]['adata'].X.sum(1).A1
        D[s]['adata'].var_names = gene_names
        D[s]['adata'].obs['unique_cell_id'] = cell_bc_seqs
        D[s]['adata'].obs['cell_names'] = cell_names
        D[s]['adata'].obs['library_id'] = np.tile(s, [D[s]['adata'].n_obs, 1])
        D[s]['adata'].uns['library_id'] = s

    return D


def load_celldata(adata, csv_filename, filter_nomatch=False):
    '''
    Adds cell annotations to the 'obs' dataframe of a ScanPy AnnData object (adata) from an imported CSV file.  
    Uses a set of unique cell identifiers (e.g. inDrops cell barcode sequences) to match cells.  These 
    identifiers are present in AnnData (in adata.obs.unique_cell_id) and in the first column of the CSV file.

    The structure of the CSV file is as follows:
    Column 1: unique cell identifiers (exact string matches to elements of adata.obs.unique_cell_id)
    Column 2: first cell annotation
    Column 3: second cell annotation
      ...          ....   
    Column n: last cell annotation  
    Column headers in the CSV file (required) will become headers of new columns in adata.obs       

    Unique cell ids in adata that no not appear in the CSV file will be annotated as 'no match'.
    'filter_nomatch' gives an option to remove these cells in the outputted version of adata.
    '''

    uID_query = adata.obs.unique_cell_id

    # load CSV header, get the names and number of IDs
    header = pd.read_csv(csv_filename, nrows=0)
    annotation_names = list(header.columns.values)[
        1:]  # ignore the first column header
    nAnnotations = len(annotation_names)

    # make a dictionary of unique cell IDs and annotations from the CSV file
    loadtxt = np.loadtxt(csv_filename, dtype='str', delimiter=',', skiprows=1)
    annotation_dict = {}
    for uID, *annots in loadtxt:   # column1 = uID, all remaining columns are annotations
        annotation_dict[uID] = annots

    # lookup each query in the dictionary, return matching annotations (or NaNs)
    annotations = []
    for j, uID in enumerate(uID_query):
        if uID in annotation_dict:
            match = annotation_dict.get(uID)
            annotations.append(match)
        else:
            annotations.append(np.repeat('no match', nAnnotations).tolist())

    # convert from list of lists to array
    annotations = np.array(annotations)

    # now copy the matched annotations to adata
    for j in range(0, nAnnotations):
        adata.obs[annotation_names[j]] = annotations[:, j]

    # if invoked, remove cells that were not present in the annotation CSV file
    if filter_nomatch:
        adata = adata[adata.obs[annotation_names[j]] != 'no match', :]

    return adata


# DATA PRE-PROCESSING

def filter_abundant_barcodes(adata, filter_cells=True, save_path='./figures/'):
    '''
    Plots a weighted histogram of transcripts per cell barcode for guiding the
    placement of a filtering threshold. Returns a filtered version of adata.  
    '''

    # If necessary, create the output directory
    if not os.path.isdir(save_path):
        os.makedirs(save_path)

    # Load counts data etc from adata
    counts = adata.obs['n_counts'].values
    threshold = adata.uns['counts_thresh']
    library_name = adata.uns['library_id']

    # Plot and format a weighted counts histogram
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(counts, bins=np.logspace(0, 6, 100), weights=counts / sum(counts))
    ax.set_xscale('log')
    ax.set_xlabel('Transcripts per cell barcode')
    ax.set_ylabel('Fraction of total transcripts')
    ax.set_title(library_name + ' (Weighted)')

    # Overlay the counts threshold as a vertical line
    ax.plot([threshold, threshold], ax.get_ylim())

    # Save figure to file
    fig.tight_layout()
    plt.savefig(save_path + 'barcode_hist_' + library_name + '.png')
    plt.show()
    plt.close()

    # Print the number of cell barcodes that will be retained vs. the total number of
    # cell barcodes in the library
    ix = counts >= threshold
    print('Filtering barcodes for', library_name,
          ' (', np.sum(ix), '/', counts.shape[0], ')')

    # Return a filtered version of adata
    if filter_cells:
        sc.pp.filter_cells(adata, min_counts=threshold, inplace=True)

    return adata, fig, ax


# VARIABLE GENES

def get_vscores(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1):
    '''
    Calculate v-score (above-Poisson noise statistic) for genes in the input counts matrix
    Return v-scores and other stats
    '''

    ncell = E.shape[0]

    mu_gene = E.mean(axis=0).A.squeeze()
    gene_ix = np.nonzero(mu_gene > min_mean)[0]
    mu_gene = mu_gene[gene_ix]

    tmp = E[:, gene_ix]
    tmp.data **= 2
    var_gene = tmp.mean(axis=0).A.squeeze() - mu_gene ** 2
    del tmp
    FF_gene = var_gene / mu_gene

    data_x = np.log(mu_gene)
    data_y = np.log(FF_gene / mu_gene)

    x, y = runningquantile(data_x, data_y, fit_percentile, nBins)
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    def gLog(input): return np.log(input[1] * np.exp(-input[0]) + input[2])
    h, b = np.histogram(np.log(FF_gene[mu_gene > 0]), bins=200)
    b = b[:-1] + np.diff(b) / 2
    max_ix = np.argmax(h)
    c = np.max((np.exp(b[max_ix]), 1))

    def errFun(b2): return np.sum(abs(gLog([x, c, b2]) - y) ** error_wt)
    b0 = 0.1
    b = scipy.optimize.fmin(func=errFun, x0=[b0], disp=False)
    a = c / (1 + b) - 1

    v_scores = FF_gene / ((1 + a) * (1 + b) + b * mu_gene)
    CV_eff = np.sqrt((1 + a) * (1 + b) - 1)
    CV_input = np.sqrt(b)

    return v_scores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b


def filter_variable_genes(E, base_ix=[], min_vscore_pctl=85, min_counts=3, min_cells=3, show_vscore_plot=False, sample_name=''):
    ''' 
    Filter genes by expression level and variability
    Return list of filtered gene indices
    '''

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = get_vscores(
        E[base_ix, :])
    ix2 = Vscores > 0
    Vscores = Vscores[ix2]
    gene_ix = gene_ix[ix2]
    mu_gene = mu_gene[ix2]
    FF_gene = FF_gene[ix2]
    min_vscore = np.percentile(Vscores, min_vscore_pctl)
    ix = (((E[:, gene_ix] >= min_counts).sum(0).A.squeeze()
           >= min_cells) & (Vscores >= min_vscore))

    if show_vscore_plot:
        import matplotlib.pyplot as plt
        x_min = 0.5 * np.min(mu_gene)
        x_max = 2 * np.max(mu_gene)
        xTh = x_min * np.exp(np.log(x_max / x_min) * np.linspace(0, 1, 100))
        yTh = (1 + a) * (1 + b) + b * xTh
        plt.figure(figsize=(8, 6))
        plt.scatter(np.log10(mu_gene), np.log10(FF_gene),
                    c=[.8, .8, .8], alpha=0.3, edgecolors='')
        plt.scatter(np.log10(mu_gene)[ix], np.log10(FF_gene)[
                    ix], c=[0, 0, 0], alpha=0.3, edgecolors='')
        plt.plot(np.log10(xTh), np.log10(yTh))
        plt.title(sample_name)
        plt.xlabel('log10(mean)')
        plt.ylabel('log10(Fano factor)')
        plt.show()

    return gene_ix[ix]


# GEPHI IMPORT & EXPORT

def export_to_graphml(adata, filename='test.graphml', directed=None):    
    import igraph as ig

    adjacency = adata.uns['neighbors']['connectivities']

    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shap[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warn('The constructed graph has only {} nodes. '
                  'Your adjacency matrix contained redundant nodes.'
                  .format(g.vcount()))
    g.write_graphml(filename)


def import_pajek_xy(adata, filename='test.net'):
    
    # first determine the number of graph nodes in *.net file
    with open(filename,'r') as file:
        nNodes = 0
        for ln,line in enumerate(file):
            if line.startswith("*Edges"):
                nNodes = ln-1

    # extract xy coordinates from *.net file
    with open(filename,'r') as file:
        lines=file.readlines()[1:nNodes+1] 
        xy = np.empty((nNodes,2))
        for ln,line in enumerate(lines):
            xy[ln,0]=(float(line.split(' ')[2]))
            xy[ln,1]=(float(line.split(' ')[3]))

    # generate ForceAtlas2 data structures and update coordinates
    sc.tl.draw_graph(adata, layout='fa', iterations=1)
    adata.obsm['X_draw_graph_fa']=xy

    return adata


# CLASSIFICATION

def train_classifiers(X, labels, PCs, gene_ind):
    '''
    Trains a series of machine learning classifiers to associate individual cells with class labels.
    Does so in a low-dimensional PCA representation of the data (PCs) over pre-defined genes (gene_ind).
    '''

    # Import sklearn classifier packages
    from sklearn.model_selection import train_test_split
    from sklearn.neural_network import MLPClassifier
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.svm import SVC
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.naive_bayes import GaussianNB
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis    


    # Subset by gene indices; project X into PCA subspace
    X_ind = X[:,gene_ind]
    PCs_ind = PCs[gene_ind,:]
    X_PCA = np.matmul(X_ind,PCs_ind)
    
    # Specify classifiers and their settings 
    classifier_names = ['NearestNeighbors', 'RandomForest', 'NeuralNet', 'LDA']
    classifiers = [KNeighborsClassifier(20, weights='distance', metric='correlation'),
                   RandomForestClassifier(n_estimators=200, random_state=802),
                   MLPClassifier(random_state=802),
                   LinearDiscriminantAnalysis()]
    
    # Split data into training and test subsets
    X_train, X_test, labels_train, labels_test = train_test_split(X_PCA, labels, test_size=0.5, random_state=802)
        
    # Build a dictionary of classifiers
    scores = []
    ClassifierDict={}
    for n,name in enumerate(classifier_names):
        clf_test = classifiers[n].fit(X_train, labels_train)
        score = clf_test.score(X_test, labels_test)
        scores.append(score)
        print(name,round(score,3))
        ClassifierDict[name]=classifiers[n].fit(X_PCA, labels)
    
    # Export classifier dictionary and subspace projection objects

    return {'Classes' : np.unique(labels),
            'Classifiers' : ClassifierDict,
    		'Classifier_Scores' : dict(zip(classifier_names, scores)), 
            'PC_Loadings' : PCs,
            'Gene_Ind' : gene_ind}
   

def predict_classes(adata, Classifier):    
    '''
    '''
    X = adata.X
    X[np.isnan(X)]=0
    PCs = Classifier['PC_Loadings']
    gene_ind = Classifier['Gene_Ind']

    # First check to see if genes match between adata and Classifier 
    adata_genes = np.array(adata.var.index) 
    classifier_genes = np.array(gene_ind.index)
    if len(classifier_genes)==len(adata_genes):
        if (classifier_genes==adata_genes).all():
            # Subset by gene indices; project X into PCA subspace
            X_ind = X[:,gene_ind]
            PCs_ind = PCs[gene_ind,:]
            X_PCA = np.matmul(X_ind,PCs_ind)
    
    else:
        # Match highly variable classifier genes to adata genes, correcting for case
        adata_genes = np.array([x.upper() for x in adata_genes])
        classifier_genes = np.array([x.upper() for x in np.array(classifier_genes[gene_ind])])
        # Get overlap
        gene_overlap, dataset_ind, classifier_ind = np.intersect1d(adata_genes,classifier_genes,return_indices=True)
        # Subset by gene indices; project X into PCA subspace
        PCs_ind = PCs[gene_ind,:]
        PCs_ind = PCs_ind[classifier_ind,:]
        X_ind = X[:,dataset_ind]
        X_PCA = np.matmul(X_ind,PCs_ind)

    # Predict class labels and probabilities for each cell, store results in adata
    for n,name in enumerate(Classifier['Classifiers']):
        adata.obs['pr_'+name] = Classifier['Classifiers'][name].predict(X_PCA)
        if hasattr(Classifier['Classifiers'][name], "predict_proba"): 
            adata.obsm['proba_'+name] = Classifier['Classifiers'][name].predict_proba(X_PCA)

    return adata


# CLUSTERING
    
def plot_confusion_matrix(labels_A, labels_B,
                          normalize=True,
                          title=None,
                          cmap=plt.cm.Blues,
                          overlay_values=False,
                          vmin=None,
                          vmax=None,
                          return_data=False):
    '''
    Plots a confusion matrix comparing two sets labels. 

    '''

    from sklearn.metrics import confusion_matrix
    from sklearn.utils.multiclass import unique_labels

    # Compute confusion matrix; 
    cm = confusion_matrix(labels_A, labels_B)
    non_empty_rows = cm.sum(axis=0)!=0
    non_empty_cols = cm.sum(axis=1)!=0
    cm = cm[:,non_empty_rows]
    cm = cm[non_empty_cols,:]
    cm = cm.T

    # Classes are the unique labels
    classes = np.unique(labels_A.append(labels_B))
    xaxis_labels = classes[non_empty_cols]
    yaxis_labels = classes[non_empty_rows]

    # Normalize by rows (label B)
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    # Set title, colorbar, and axis names
    if normalize:
        colorbar_label = 'Fraction Overlap'
        if not title:
            title = 'Normalized confusion matrix'
    else:
        colorbar_label = '# Overlaps'
        if not title:
        	title = 'Confusion matrix, without normalization'  
  
    if hasattr(labels_A, 'name'):
        labels_A_name = labels_A.name #.capitalize()   	
    else:
        labels_A_name = 'Label A'
    if hasattr(labels_B, 'name'):
        labels_B_name = labels_B.name #.capitalize()    	
    else:
        labels_B_name = 'Label B'

    # Generate and format figure axes
    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)

    ax.grid(False)
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           xticklabels=xaxis_labels, yticklabels=yaxis_labels,
           title=title,
           ylabel=labels_B_name,
           xlabel=labels_A_name)

    # Format tick labels
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", va='top',
             rotation_mode='anchor',fontsize=10)
    plt.setp(ax.get_yticklabels(), fontsize=10)

    # Format colorbar
    cb=ax.figure.colorbar(im, ax=ax, shrink=0.5)
    cb.ax.tick_params(labelsize=10) 
    cb.ax.set_ylabel(colorbar_label, rotation=90)
    
    # Loop over data dimensions and create text annotations
    if overlay_values:
        fmt = '.1f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, format(cm[i, j], fmt),
                        ha="center", va="center",
                        color="white" if cm[i, j] > thresh else "black",
                        size=8)
    ax.set_aspect('equal') 
    
    if return_data:
        return fig, ax, cm, xaxis_labels, yaxis_labels                        
    else:
        return fig, ax                   


# DIFFERENTIAL EXPRESSION

def get_dynamic_genes(adata, sliding_window=100, fdr_alpha = 0.05):

    # Input an AnnData object that has already been subsetted to cells and genes of interest.
    # Cells are ranked by dpt pseudotime. Genes are tested for significant differential expression 
    # between two sliding windows corresponding the highest and lowest average expression. FDR values
    # are then calculated by thresholding p-values calculated from randomized data.
    # Returns a copy of adata with the following fields added: 
    #   adata.var['dyn_peak_cell']: pseudotime-ordered cell with the highest mean expression
    #   adata.var['dyn_fdr']: fdr-corrected p-value for differential expression
    #   adata.var['dyn_fdr_flag']: boolean flag, true if fdr <= fdr_alpha

    import scipy.stats

    # Function for calculating p-values for each gene from min & max sliding window expression values
    def get_slidingwind_pv(X, sliding_window):
        # construct a series of sliding windows over the cells in X
        wind=[]
        nCells = X.shape[0]
        for k in range(nCells-sliding_window+1):    
            wind.append(list(range(k, k+sliding_window)))
        # calculate p-values on the sliding windows
        pv = []
        max_cell_this_gene = []
        nGenes = X.shape[1]
        for j in range(nGenes):
            tmp_X_avg = []
            # get mean expression of gene j in each sliding window k
            for k in range(len(wind)-1):    
                tmp_X_avg.append(np.mean(X[wind[k],j]))
            # determine min and max sliding windows for this gene
            max_wind = np.argmax(tmp_X_avg)
            min_wind = np.argmin(tmp_X_avg)
            # determine if this gene displays significant differential expression
            _,p=scipy.stats.ttest_ind(X[wind[max_wind],j],X[wind[min_wind],j])
            pv.append(p[0])
            max_cell_this_gene.append(max_wind)
        return np.array(pv), np.array(max_cell_this_gene)

    # import counts and pseudotime from the AnnData object
    nCells = adata.shape[0]
    nGenes = adata.shape[1]
    cell_order = np.argsort(adata.obs['dpt_pseudotime'])
    if scipy.sparse.issparse(adata.X):
        X = adata.X[cell_order,:].todense()
    else:
        X = adata.X[cell_order,:]

    # calculate p values on the pseudotime-ordered data
    pv, peak_cell = get_slidingwind_pv(X, sliding_window)
    adata.var['dyn_peak_cell'] = peak_cell#np.argsort(gene_ord)
    print('done calculating p-values')
    
    # calculate p values on the randomized data
    np.random.seed(802)
    X_rand = X[np.random.permutation(cell_order),:]
    pv_rand, _ = get_slidingwind_pv(X_rand, sliding_window)
    print('done calculating randomized p-values')

    # calculate fdr as the fraction of randomized p-values that exceed this p-value
    fdr = []
    fdr_flag = []
    for j in range(nGenes):
        fdr.append(sum(pv_rand <= pv[j])/nGenes)
        fdr_flag.append(fdr[j] <= fdr_alpha)
    adata.var['dyn_fdr'] = fdr
    adata.var['dyn_fdr_flag'] = fdr_flag
    print('done calculating fdr')

    return adata

    
# PLOTTING

def format_axes(eq_aspect='all', rm_colorbar=False):
    '''
    Gets axes from the current figure and applies custom formatting options
    In general, each parameter is a list of axis indices (e.g. [0,1,2]) that will be modified
    Colorbar is assumed to be the last set of axes
    '''
    
    # get axes from current figure
    ax = plt.gcf().axes

    # format axes aspect ratio
    if eq_aspect is not 'all':
        for j in eq_aspect:
            ax[j].set_aspect('equal') 
    else:
        for j in range(len(ax)):
            ax[j].set_aspect('equal') 

    # remove colorbar
    if rm_colorbar:
        j=len(ax)-1
        if j>0:
            ax[j].remove()

# SCANPY 
# score_genes function with random_state behavior fixed

def score_genes(
        adata,
        gene_list,
        ctrl_size=50,
        gene_pool=None,
        n_bins=25,
        score_name='score',
        random_state=0,
        copy=False,
        use_raw=False):  
    adata = adata.copy() if copy else adata

    np.random.seed(random_state)

    gene_list_in_var = []
    var_names = adata.raw.var_names if use_raw else adata.var_names
    for gene in gene_list:
        if gene in var_names:
            gene_list_in_var.append(gene)
    gene_list = set(gene_list_in_var[:])

    if not gene_pool:
        gene_pool = list(var_names)
    else:
        gene_pool = [x for x in gene_pool if x in var_names]

    _adata = adata.raw if use_raw else adata
    if scipy.sparse.issparse(_adata.X):
        obs_avg = pd.Series(
            np.nanmean(
                _adata[:, gene_pool].X.toarray(), axis=0), index=gene_pool)  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(_adata[:, gene_pool].X, axis=0), index=gene_pool)  # average expression of genes

    obs_avg = obs_avg[np.isfinite(obs_avg)] # Sometimes (and I don't know how) missing data may be there, with nansfor

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    control_genes = set()

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))  # uses full r_genes if ctrl_size > len(r_genes)

    # To index, we need a list - indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)


    X_list = _adata[:, gene_list].X
    if scipy.sparse.issparse(X_list): X_list = X_list.toarray()
    X_control = _adata[:, control_genes].X
    if scipy.sparse.issparse(X_control): X_control = X_control.toarray()
    X_control = np.nanmean(X_control, axis=1)

    if len(gene_list) == 0:
        return adata if copy else None
    elif len(gene_list) == 1:
        score = _adata[:, gene_list].X - X_control
    else:
        score = np.nanmean(X_list, axis=1) - X_control

    adata.obs[score_name] = pd.Series(np.array(score).ravel(), index=adata.obs_names)

    return adata if copy else None

