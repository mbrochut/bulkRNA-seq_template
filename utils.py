import pandas as pd
import plotly.express as px
import numpy as np
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import seaborn as sns
import textalloc as ta


from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list

import plotly.graph_objects as go
from plotly.colors import DEFAULT_PLOTLY_COLORS


def plot_mean_over_variance(adata,save=None):
    # Calculate mean and variance for each gene (row-wise)
    mean_counts = np.mean(adata.X, axis=0)
    variance_counts = np.var(adata.X, axis=0)

    # Create a dataframe to store mean and variance values
    df = pd.DataFrame({'mean_counts': mean_counts, 'variance_counts': variance_counts})

    # Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df['mean_counts'], df['variance_counts'], alpha=0.6)
    plt.plot(df['mean_counts'], df['mean_counts'], color='red', label='y = x')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Mean Counts (log scale)')
    plt.ylabel('Variance Counts (log scale)')
    plt.title('Variance over Mean Plot')
    plt.legend()
    plt.grid(True, which="both", ls="--", lw=0.5)
    if save:
        plt.savefig(f"./figures/QC/{save}png")
    plt.show()

def plot_sample_distance_heatmap(
    adata,
    layer: str | None = None,
    metric: str = "euclidean",
    method: str = "ward",
    figsize: tuple = (10, 10),
    cmap: str = "rocket",
    save: str | None = None,
):
    """
    Plot a clustered heatmap of sample-to-sample distances.

    Parameters
    ----------
    adata : AnnData
        AnnData object with samples in observations.
    layer : str or None, default None
        Layer to use. If None, use adata.X.
    metric : str, default "euclidean"
        Distance metric for pdist.
    method : str, default "ward"
        Linkage method for hierarchical clustering.
    figsize : tuple, default (10, 10)
        Figure size.
    cmap : str, default "rocket"
        Seaborn colormap.
    save : str or None, default None
        Filename (without extension) to save the plot in ./results/QC/.
    """

    # Select matrix
    X = adata.layers[layer] if layer is not None else adata.X
    X = X.toarray() if hasattr(X, "toarray") else X

    # Compute distances
    sample_dists = pdist(X, metric=metric)
    dist_matrix = squareform(sample_dists)

    sample_names = adata.obs_names.to_list()
    dist_df = pd.DataFrame(dist_matrix, index=sample_names, columns=sample_names)

    # Linkage
    row_linkage = linkage(sample_dists, method=method)
    col_linkage = linkage(sample_dists, method=method)

    # Plot
    g = sns.clustermap(
        dist_df,
        row_cluster=True,
        col_cluster=True,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap=sns.color_palette(cmap, as_cmap=True),
        figsize=figsize,
    )

    if save:
        plt.savefig(f"./figures/QC/{save}.png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_raw_vs_normed_counts(
    adata,
    raw_layer: str = "raw_counts",
    normed_layer: str = "normed_counts",
    figsize: tuple = (15, 10),
    save: str | None = None,
):
    """
    Plot boxplots of log10-transformed raw and normalized counts per sample.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing count layers.
    raw_layer : str, default "raw_counts"
        Layer name for raw counts.
    normed_layer : str, default "normed_counts"
        Layer name for normalized counts.
    figsize : tuple, default (15, 10)
        Figure size.
    save : str or None, default None
        Filename (without extension) to save the plot in ./results/QC/.
    """

    def _prepare_long_df(X, gene_names, sample_names):
        X = X.toarray() if hasattr(X, "toarray") else X
        X = np.log10(X + 1)
        df = pd.DataFrame(X.T, index=gene_names, columns=sample_names)
        return df.melt(var_name="Sample", value_name="Counts")

    raw_long = _prepare_long_df(
        adata.layers[raw_layer],
        adata.var_names,
        adata.obs_names,
    )

    normed_long = _prepare_long_df(
        adata.layers[normed_layer],
        adata.var_names,
        adata.obs_names,
    )

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=figsize, sharex=True)

    sns.boxplot(data=raw_long, x="Sample", y="Counts", ax=axes[0])
    axes[0].set_title("Raw Counts per Sample")
    axes[0].tick_params(axis="x", rotation=90)

    sns.boxplot(data=normed_long, x="Sample", y="Counts", ax=axes[1], color="red")
    axes[1].set_title("Normalized Counts per Sample")
    axes[1].tick_params(axis="x", rotation=90)

    plt.tight_layout()

    if save:
        plt.savefig(f"./figures/QC/{save}.png", dpi=300, bbox_inches="tight")

    plt.show()


def volcano_plotly(df,
    cutoff_fc=np.log2(1.5),
    cutoff_p=0.05,
    highlight_n=10,
    height=900,
    width=1000,
    col_gene_name = 'gene_name',
    color_palette=None,
    xlim=(-8,8)
    ):
    """
    Function to create a volcano plot with Plotly and label the most important genes.
    
    Parameters:
    df : DataFrame
        DataFrame containing columns 'log2FoldChange', 'padj', and 'gene_name'.
    cutoff_fc : float, optional
        Fold-change cutoff for determining significance. Default is 1.
    cutoff_p : float, optional
        p-value cutoff for determining significance. Default is 0.05.
    highlight_n : int, optional
        Number of genes to highlight on the plot by labeling the most important genes.
    height : int, optional
        Height of the plot. Default is 900.
    width : int, optional
        Width of the plot. Default is 1000.
        
    Returns:
    Plotly figure object
    """
    if color_palette:
        color_palette = color_palette
    else:
        color_palette = ('#00BFC4','#F8766D')

    # Create new columns for colors and labels based on the cutoff criteria
    df['keyvals.colour'] = np.where(
        (df['log2FoldChange'] < -cutoff_fc) & (df['padj'] < cutoff_p), color_palette[1], 
        np.where((df['log2FoldChange'] > cutoff_fc) & (df['padj'] < cutoff_p), color_palette[0], 'grey')
    )

    df['keyvals.label'] = np.where(
        (df['log2FoldChange'] < -cutoff_fc) & (df['padj'] < cutoff_p), 'Down regulated', 
        np.where((df['log2FoldChange'] > cutoff_fc) & (df['padj'] < cutoff_p), 'Up regulated', 'NS')
    )

    # Add a new column for -log10 of the padj values
    df['-log10(padj)'] = -np.log10(df['padj'])

    # Identify the most important genes (e.g., based on top/bottom n log2FoldChange)
    most_important_genes = df[(df['padj'] < cutoff_p) & (df['log2FoldChange'].abs() >cutoff_fc)].sort_values(by='padj').head(highlight_n)

    # Create the scatter plot
    trace = go.Scatter(
        x=df['log2FoldChange'],
        y=df['-log10(padj)'],
        mode='markers',
        marker=dict(color=df['keyvals.colour']),
        text=df.apply(lambda row: f"Gene: {row[col_gene_name]}<br>log2FoldChange: {row['log2FoldChange']}<br>padj: {row['padj']}", axis=1),
        hoverinfo='text'
    )

    # Add text labels for the most important genes, while still having hover info
    trace_labels = go.Scatter(
        x=most_important_genes['log2FoldChange'],
        y=most_important_genes['-log10(padj)'],
        mode='markers+text',
        text=most_important_genes[col_gene_name],
        textposition='top center',
        marker=dict(color=most_important_genes['keyvals.colour']),  # Larger markers for labeled points
        hoverinfo='text',
        textfont=dict(size=12),  # Set font size for gene labels
        texttemplate='%{text}',  # Ensure only gene names are shown as static labels
        customdata=most_important_genes.apply(lambda row: f"Gene: {row[col_gene_name]}<br>log2FoldChange: {row['log2FoldChange']}<br>padj: {row['padj']}", axis=1),
        hovertemplate='%{customdata}'  # Custom hover template for labeled points
    )

    # Define the maximum Y-value for the plot (for cutoff line lengths)
    maxValue = df['-log10(padj)'].max()

    # Create cutoff lines
    cutoff_lines = [
        # Vertical lines for fold-change cutoff
        dict(type="line", x0=cutoff_fc, x1=cutoff_fc, y0=0, y1=maxValue, line=dict(dash="dash", color="black")),
        dict(type="line", x0=-cutoff_fc, x1=-cutoff_fc, y0=0, y1=maxValue, line=dict(dash="dash", color="black")),
        
        # Horizontal line for p-value cutoff
        dict(type="line", x0=0, x1=1, xref='paper', y0=-np.log10(cutoff_p), y1=-np.log10(cutoff_p), line=dict(dash="dash", color="black"))
    ]

    # Create the layout
    layout = go.Layout(
        title="Volcano Plot",
        xaxis=dict(title='log2FoldChange'),
        yaxis=dict(title='-log10(padj)'),
        height=height,
        width=width,
        showlegend=False,
        shapes=cutoff_lines,
    )

    # Create the figure
    fig = go.Figure(data=[trace, trace_labels], layout=layout)

    xmin = max(fig.data[0].x.min(), xlim[0])
    xmax = min(fig.data[0].x.max(), xlim[1])
 
    fig.update_xaxes(range=[xmin, xmax])
    # Display the plot
    fig.show(config={"responsive": True})

def volcano_matplotlib_with_textalloc(
    df,
    cutoff_fc=0.58,
    cutoff_p=0.05,
    highlight_n=10,
    height=9,
    width=10,
    x_range=None,
    color_palette=None,
    ):
    """
    Function to create a volcano plot with Matplotlib and label the most important genes using textalloc.

    Parameters:
    df : DataFrame
        DataFrame containing columns 'log2FoldChange', 'padj', and 'gene_name'.
    cutoff_fc : float, optional
        Fold-change cutoff for determining significance. Default is 1.
    cutoff_p : float, optional
        p-value cutoff for determining significance. Default is 0.001.
    highlight_n : int, optional
        Number of genes to highlight on the plot by labeling the most important genes.
    height : float, optional
        Height of the plot in inches. Default is 9.
    width : float, optional
        Width of the plot in inches. Default is 10.
    x_range : tuple, optional
        Range for the x-axis (log2FoldChange). Default is None (auto).

    Returns:
    Matplotlib figure and axis
    """
    if color_palette:
        color_palette = color_palette
    else:
        color_palette = ('#00BFC4','#F8766D')

    # Prepare data for coloring and labeling
    df['keyvals.colour'] = np.where(
        (df['log2FoldChange'] < -cutoff_fc) & (df['padj'] < cutoff_p), color_palette[1],
        np.where((df['log2FoldChange'] > cutoff_fc) & (df['padj'] < cutoff_p), color_palette[0], 'grey')
    )
    df['keyvals.label'] = np.where(
        (df['log2FoldChange'] < -cutoff_fc) & (df['padj'] < cutoff_p), 'Down regulated',
        np.where((df['log2FoldChange'] > cutoff_fc) & (df['padj'] < cutoff_p), 'Up regulated', 'NS')
    )
    df['-log10(padj)'] = -np.log10(df['padj'])

    # Identify the most important genes
    most_important_genes = df[(df['padj'] < cutoff_p) & (df['log2FoldChange'].abs() > cutoff_fc)]\
        .sort_values(by='padj').head(highlight_n)

    # Initialize the plot
    plt.figure(figsize=(width, height))
    sns.set(style="white")

    # Scatter plot
    plt.scatter(df['log2FoldChange'], df['-log10(padj)'], c=df['keyvals.colour'], alpha=0.7, edgecolor=None)

    # Highlight most important genes
    plt.scatter(most_important_genes['log2FoldChange'], most_important_genes['-log10(padj)'],
                c=most_important_genes['keyvals.colour'], edgecolor='black', s=80, label='Highlighted')

    # Add cutoff lines
    plt.axhline(-np.log10(cutoff_p), color='black', linestyle='--', linewidth=0.8, label=f'Cutoff: p={cutoff_p}')
    plt.axvline(cutoff_fc, color='black', linestyle='--', linewidth=0.8, label=f'Cutoff: FC={cutoff_fc}')
    plt.axvline(-cutoff_fc, color='black', linestyle='--', linewidth=0.8)

    # Use textalloc for better label placement
    ta.allocate(
        plt.gca(),
        most_important_genes['log2FoldChange'].values,
        most_important_genes['-log10(padj)'].values,
        most_important_genes['gene_name'].values,
        x_scatter=df['log2FoldChange'].values,
        y_scatter=df['-log10(padj)'].values,
        max_distance=0.07,
        #avoid_label_lines_overlap = True,
        #direction = 'north',
    )

    # Customize plot
    plt.title("Volcano Plot", fontsize=28)
    plt.xlabel('log2FoldChange', fontsize=22)
    plt.ylabel('-log10(padj)', fontsize=22)
    if x_range:
        plt.xlim(x_range)
    plt.ylim(0, df['-log10(padj)'].max() + 1)
    plt.legend(loc='upper left', fontsize=14)
    plt.tight_layout()

    # Show plot
    plt.show()


def ORA_plotly(ora_res,width=1000,height=1000,topN=30):
    to_display = ora_res.head(topN)
    to_display['-log10pvalue'] = -np.log10(to_display['P-value'])
    to_display = to_display.sort_values(by='Count',ascending=True)

    fig = px.scatter(
    to_display,
    x='Count',
    y='Term',
    size='Overlap_Ratio',
    color='P-value',
    color_continuous_scale='Viridis',
    title='ORA Dot Plot',
)

    # Customize the plot
    fig.update_layout(
        xaxis_title='Count',
        yaxis_title='Pathway',
        coloraxis_colorbar=dict(
        title='P-value',
        tickformat=".1e"  # Use scientific notation for colorbar ticks
        ),
        showlegend=False,
        width=width,
        height=height,
    )
    # Customize the size legend (dot size legend)
   
    # Display the plot
    fig.show()


def GSEA_plotly(gsea_res, width=1000, height=1000):
    # Prepare the data to display, selecting the top 30 entries and sorting by '-log10Qvalue'
    to_display = gsea_res.head(30)
    to_display = to_display.sort_values(by='-log10Qvalue', ascending=True)

    # Determine if there are negative NES values
    if min(to_display['NES']) <0 and max(to_display['NES']) >0:
        palette = px.colors.diverging.RdBu[::-1]
        color_range  = [min(to_display['NES']),max(to_display['NES'])]
    elif min(to_display['NES']) <0 and max(to_display['NES']) <0:
        palette = px.colors.sequential.Blues[::-1]
        color_range=[min(to_display['NES']),0]
    else:
        palette = px.colors.sequential.Reds
        color_range=[0,max(to_display['NES'])]

    # Create the scatter plot using Plotly Express
    fig = px.scatter(
        to_display,
        x='-log10Qvalue',
        y='Term',
        size='Tag_%',  # Dot size determined by the 'Tag_%' column
        color='NES',   # Color determined by the 'NES' column
        color_continuous_scale=palette,  # Inverted RdBu color scale
        range_color=color_range,  # Dynamically set the color scale range
        title='GSEA Dot Plot',
    )

    # Customize the plot layout
    fig.update_layout(
        xaxis_title='-log10Qvalue',
        yaxis_title='Pathway',
        coloraxis_colorbar=dict(title='NES'),  # Change colorbar title to 'NES'
        showlegend=False,
        width=width,
        height=height,
    )

    # Set x-axis to start from 0
    fig.update_xaxes(range=[0, None])

    # Add a red dashed vertical line at -log10(0.05) ≈ 1.3
    significance_threshold = -np.log10(0.05)
    fig.add_shape(
        type='line',
        x0=significance_threshold, x1=significance_threshold,
        y0=0, y1=1,
        xref='x', yref='paper',  # yref='paper' makes the line span the entire y-axis
        line=dict(color='red', width=2, dash='dash')
    )

    # Display the plot
    fig.show()

def confidence_ellipse(x, y, n_std=1.96, size=100):
    """
    Get the covariance confidence ellipse of *x* and *y*.
    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.
    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.
    size : int
        Number of points defining the ellipse
    Returns
    -------
    String containing an SVG path for the ellipse
    
    References (H/T)
    ----------------
    https://matplotlib.org/3.1.1/gallery/statistics/confidence_ellipse.html
    https://community.plotly.com/t/arc-shape-with-path/7205/5
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    theta = np.linspace(0, 2 * np.pi, size)
    ellipse_coords = np.column_stack([ell_radius_x * np.cos(theta), ell_radius_y * np.sin(theta)])
    
    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    x_scale = np.sqrt(cov[0, 0]) * n_std
    x_mean = np.mean(x)

    # calculating the stdandard deviation of y ...
    y_scale = np.sqrt(cov[1, 1]) * n_std
    y_mean = np.mean(y)
  
    translation_matrix = np.tile([x_mean, y_mean], (ellipse_coords.shape[0], 1))
    rotation_matrix = np.array([[np.cos(np.pi / 4), np.sin(np.pi / 4)],
                                [-np.sin(np.pi / 4), np.cos(np.pi / 4)]])
    scale_matrix = np.array([[x_scale, 0],
                            [0, y_scale]])
    ellipse_coords = ellipse_coords.dot(rotation_matrix).dot(scale_matrix) + translation_matrix
        
    path = f'M {ellipse_coords[0, 0]}, {ellipse_coords[0, 1]}'
    for k in range(1, len(ellipse_coords)):
        path += f'L{ellipse_coords[k, 0]}, {ellipse_coords[k, 1]}'
    path += ' Z'
    return path


def pca_from_adata(adata):
    # PCA Results
    explained_variance = adata.uns['pca']['variance_ratio'][:2]  # Variance explained by PCs
    adata.obs[['PC1','PC2']] =  adata.obsm['X_pca'][:,:2]
    df_pca = adata.obs
    return df_pca, explained_variance


def confidence_ellipse(x, y, n_std=1.96, size=100):
    """
    Get the covariance confidence ellipse of *x* and *y*.
    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.
    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.
    size : int
        Number of points defining the ellipse
    Returns
    -------
    String containing an SVG path for the ellipse
    
    References (H/T)
    ----------------
    https://matplotlib.org/3.1.1/gallery/statistics/confidence_ellipse.html
    https://community.plotly.com/t/arc-shape-with-path/7205/5
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    theta = np.linspace(0, 2 * np.pi, size)
    ellipse_coords = np.column_stack([ell_radius_x * np.cos(theta), ell_radius_y * np.sin(theta)])
    
    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    x_scale = np.sqrt(cov[0, 0]) * n_std
    x_mean = np.mean(x)

    # calculating the stdandard deviation of y ...
    y_scale = np.sqrt(cov[1, 1]) * n_std
    y_mean = np.mean(y)
  
    translation_matrix = np.tile([x_mean, y_mean], (ellipse_coords.shape[0], 1))
    rotation_matrix = np.array([[np.cos(np.pi / 4), np.sin(np.pi / 4)],
                                [-np.sin(np.pi / 4), np.cos(np.pi / 4)]])
    scale_matrix = np.array([[x_scale, 0],
                            [0, y_scale]])
    ellipse_coords = ellipse_coords.dot(rotation_matrix).dot(scale_matrix) + translation_matrix
        
    path = f'M {ellipse_coords[0, 0]}, {ellipse_coords[0, 1]}'
    for k in range(1, len(ellipse_coords)):
        path += f'L{ellipse_coords[k, 0]}, {ellipse_coords[k, 1]}'
    path += ' Z'
    return path


def plot_pca_plotly(df_pca, explained_variance, title='',save=None):
    """
    Create a Plotly PCA plot with custom styling.

    Parameters:
    df_pca (pd.DataFrame): DataFrame containing PCA data with columns 'PC1', 'PC2', and 'condition'.
    explained_variance (list): List containing the explained variance for each principal component.
    fig_to_save (str): Optional filename to save the figure.

    Returns:
    None
    """
    # Plotly Figure
    fig = go.Figure()
    custome_order = df_pca['condition'].unique().tolist()
    custome_order = [custome_order[1],custome_order[0]]+custome_order[2:]
    for target_value, target_name in enumerate(custome_order):
        color = DEFAULT_PLOTLY_COLORS[target_value]
        data = df_pca.loc[df_pca['condition'] == target_name, :]

        fig.add_trace(
            go.Scatter(
                x=data['PC1'],
                y=data['PC2'],
                name=target_name,
                mode='markers',  # Remove text mode
                marker={'color': color, 'size': 25},
            )
        )

        fig.add_shape(
            type='path',
            path=confidence_ellipse(data['PC1'], data['PC2']),
            line={'dash': 'dot'},
            line_color=color
        )

    fig.update_layout(
        title=f"{title}",
        legend_title="Condition",
        width=900,  # Adjust width here
        height=800,  # Adjust height here
        paper_bgcolor='white',  # Outer background color
        plot_bgcolor='white',   # Inner plot area background color
        xaxis=dict(
            mirror=True,
            title=f"PC1 ({round(explained_variance[0] * 100)}%)",
            showgrid=False,       # Hide grid lines
            showline=True,
            title_font=dict(size=18),
            linecolor='black',
            linewidth=2           # Thick border
        ),
        yaxis=dict(
            mirror=True,
            title=f"PC2 ({round(explained_variance[1] * 100)}%)",
            showgrid=False,       # Hide grid lines
            showline=True,
            title_font=dict(size=18),
            linecolor='black',
            linewidth=2           # Thick border
        ),
        font=dict(size=18)
    )

    if save:
        fig.write_image(f"./figures/QC/{save}.png")
    fig.show()


def prepare_df_for_top_N_pathways(df,db = 'GO',db_col='db',p_value_col='FDR q-val',contrast_col='contrast',top_n = 10):
    
    full_df_DB = df[df[db_col]==db]
    top_n_per_contrast = full_df_DB.sort_values([p_value_col]).groupby([contrast_col]).head(top_n)
    top_n_per_contrast = top_n_per_contrast.drop_duplicates('Term')
    pathways_of_interest = top_n_per_contrast['Term'].tolist()
    db_ready_to_plot = df[df['Term'].isin(pathways_of_interest)]
    return db_ready_to_plot



 # Define a function to add significance stars based on the p-value
def significance_stars(pval):
    if pval < 0.0001:
        return '****'
    elif pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return ''


def plot_NES_heatmap_with_significance(
    df,
    custom_order=None,
    width=1000,
    height=1500
):
    # --------------------------------------------------
    # 0. Edge cases
    # --------------------------------------------------
    if df is None or df.empty:
        print("Empty dataframe: nothing to plot.")
        return None

    df = df.copy()

    # --------------------------------------------------
    # 1. Ordering
    # --------------------------------------------------
    if custom_order is None:
        custom_order = df['contrast'].unique().tolist()

    df['contrast'] = pd.Categorical(
        df['contrast'],
        categories=custom_order,
        ordered=True
    )
    df = df.sort_values(by=['contrast', 'Term'])

    # --------------------------------------------------
    # 2. Pivot
    # --------------------------------------------------
    nes_matrix = df.pivot_table(index='Term', columns='contrast', values='NES')
    fdr_matrix = df.pivot_table(index='Term', columns='contrast', values='FDR q-val non zero')

    if nes_matrix.empty or fdr_matrix.empty:
        print("Empty pivot table: nothing to plot.")
        return None

    # --------------------------------------------------
    # 3. Transform for clustering
    # --------------------------------------------------
    log_fdr_matrix = -np.log10(fdr_matrix)
    log_fdr_matrix.replace([np.inf, -np.inf], np.nan, inplace=True)

    # --------------------------------------------------
    # 4. Clustering (only if >1 row AND >1 column)
    # --------------------------------------------------
    if log_fdr_matrix.shape[0] > 1 and log_fdr_matrix.shape[1] > 1:
        Z = linkage(log_fdr_matrix.fillna(0), method='ward')
        clustered_terms = leaves_list(Z)

        nes_matrix_clustered = nes_matrix.iloc[clustered_terms]
        fdr_matrix_clustered = fdr_matrix.iloc[clustered_terms]
    else:
        nes_matrix_clustered = nes_matrix
        fdr_matrix_clustered = fdr_matrix

    # --------------------------------------------------
    # 5. Significance stars
    # --------------------------------------------------
    stars_matrix = fdr_matrix_clustered.applymap(significance_stars)

    # --------------------------------------------------
    # 6. Plot
    # --------------------------------------------------
    fig = px.imshow(
        nes_matrix_clustered,
        color_continuous_scale='RdBu_r',
        aspect='auto',
        labels=dict(x="Condition", y="Pathway", color="NES"),
        title="Top pathways per contrast (GSEA)",
        zmin=nes_matrix.min().min(),
        zmax=nes_matrix.max().max(),
        text_auto=False,  # we control text manually
    )

    # Add stars
    fig.data[0]['text'] = stars_matrix.values
    fig.data[0]['texttemplate'] = '%{text}'

    fig.update_layout(
        xaxis_title="Condition",
        yaxis_title="Pathway",
        width=width,
        height=height,
    )

    fig.show()


def plot_pathway_enrichment_heatmap(
    df,
    contrast_col='contrast',
    p_value='Adjusted P-value',
    custom_order=None,
    width=1000,
    height=1500,
    zmin=None,
    zmax=None
):
    # --------------------------------------------------
    # 0. Edge cases
    # --------------------------------------------------
    if df is None or df.empty:
        print("Empty dataframe: nothing to plot.")
        return None

    # --------------------------------------------------
    # 1. Ordering
    # --------------------------------------------------
    if custom_order is None:
        custom_order = df[contrast_col].unique().tolist()

    df = df.copy()
    df[contrast_col] = pd.Categorical(
        df[contrast_col],
        categories=custom_order,
        ordered=True
    )
    df = df.sort_values(by=[contrast_col, 'Term'])

    # --------------------------------------------------
    # 2. Pivot
    # --------------------------------------------------
    heatmap_data = df.pivot_table(
        index='Term',
        columns=contrast_col,
        values=p_value
    )

    if heatmap_data.empty:
        print("Empty pivot table: nothing to plot.")
        return None

    # --------------------------------------------------
    # 3. Transform
    # --------------------------------------------------
    heatmap_data_log = -np.log10(heatmap_data)
    heatmap_data_log.replace([np.inf, -np.inf], np.nan, inplace=True)

    # --------------------------------------------------
    # 4. Clustering (only if >1 row)
    # --------------------------------------------------
    if heatmap_data_log.shape[0] > 1:
        Z = linkage(heatmap_data_log.fillna(0), method='ward')
        clustered_terms = leaves_list(Z)
        heatmap_data_clustered = heatmap_data_log.iloc[clustered_terms]
        heatmap_data_clustered_original = heatmap_data.iloc[clustered_terms]
    else:
        # No clustering
        heatmap_data_clustered = heatmap_data_log
        heatmap_data_clustered_original = heatmap_data

    # --------------------------------------------------
    # 5. Plot
    # --------------------------------------------------
    fig = px.imshow(
        heatmap_data_clustered,
        color_continuous_scale='Viridis',
        aspect='auto',
        labels=dict(
            x="Condition",
            y="Pathway",
            color=f"-log10({p_value})"
        ),
        title="Pathway Enrichment Heatmap",
        zmin=zmin,
        zmax=zmax,
    )

    # Hover info
    fig.data[0]['hovertemplate'] = (
        'Condition: %{x}<br>' +
        'Pathway: %{y}<br>' +
        '-log10(FDR q-val): %{z}<br>' +
        'Original FDR q-val: %{customdata}<extra></extra>'
    )

    fig.data[0]['customdata'] = heatmap_data_clustered_original.values

    fig.update_layout(
        xaxis_title="Condition",
        yaxis_title="Pathway",
        width=width,
        height=height,
    )

    fig.show()
