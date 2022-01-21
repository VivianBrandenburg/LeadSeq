import seaborn as sns
import matplotlib.pyplot as plt



colors_validation = {'25 °C':'#557499', '37 °C':'#ff6d6d'}
colors_l25 = {'unpaired':'#557499', 'paired':'#557499'}
colors_l37 = {'unpaired':'#ff6d6d', 'paired':'#ff6d6d'}
colors_SDS = {'candidate':'#a8c179', 'RBS':'#557499'}



def plot_validation(data, outfile):    
    sns.set(rc={'figure.figsize':(4.55,6)})
    sns.set_style('whitegrid')
    sns.set_context('talk')

    sns.barplot(data=data, y='variable', x='value', hue='sample', capsize = 0.1 , errwidth=1.2, order = ['selectivity', 'sensitivity'], palette = colors_validation)
    plt.xlim(0,1)
    plt.xlabel('')
    plt.ylabel('')
    plt.savefig(outfile, dpi=300, bbox_inches = 'tight')
    plt.show()
    
    
    

def plot_score_distribution(data, label, outfile):
    sns.set(rc={'figure.figsize':(6,6)})
    sns.set_style('whitegrid')
    sns.set_context('talk')
    
    color_palette = colors_l37 if label=='l37' else colors_l25
    
    sns.boxplot(data = data, y='scores', x='paired', palette = color_palette)
    plt.ylim(-0.1,11)
    plt.xlabel('')
    plt.ylabel('lead scores')
    plt.title('lead scores at '+ label.replace('l','')+ ' °C')
    plt.savefig(outfile, dpi=300, bbox_inches = 'tight')
    plt.show()
  


def plot_SDS(data, outfile):
    sns.set(font_scale=1.3, rc={'figure.figsize':(6,5)})
    sns.set_style('whitegrid')
    sns.set_context('talk')
    
    sns.scatterplot(data=data, x='delta', y='mwu', hue='label', 
                    hue_order = ['RBS', 'candidate'], 
                    palette=colors_SDS, s=100, alpha=0.6, linewidth=0.2)
    plt.axhline(y=0, color='grey', linestyle='-', linewidth=1)
    plt.axvline(x=0, color='grey', linestyle='-', linewidth=1)
    plt.savefig(outfile, dpi=300, bbox_inches = 'tight')
    plt.show()
    