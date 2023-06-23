import logomaker as lm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def ml_sample01():
    # create dataframe with count of each nucleotide at each position in the motif
    counts_df = pd.DataFrame({
        'A': [80, 5, 7, 49, 23, 35, 17],
        'C': [5, 7, 73, 34, 18, 25, 7],
        'G': [7, 79, 8, 7, 19, 11, 28],
        'T': [8, 9, 12, 10, 20, 9, 18]
    })

    print(counts_df)

    # convert counts to probabilities
    probs_df = counts_df.apply(lambda x: x / x.sum(), axis=1)

    print(probs_df)

    # create Logo object
    logo = lm.Logo(probs_df, color_scheme='classic',
                   width=.8, stack_order='small_on_top',
                   font_name='Arial')

    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.ax.set_title('DNA Binding Motif', fontsize=16)
    logo.ax.set_ylabel('Frequency', fontsize=14)
    logo.ax.set_xlabel('Position', fontsize=14)

    plt.show()

def ml_sample02():
    # create dataframe with count of each nucleotide at each position in the motif
    counts_df = pd.DataFrame({
        'A': [80, 5, 7, 49, 23, 35, 17],
        'C': [5, 7, 73, 34, 18, 25, 7],
        'G': [7, 79, 8, 7, 19, 11, 28],
        'T': [8, 9, 12, 10, 20, 9, 18]
    })

    print(counts_df)

    freq_matrix = counts_df.div(counts_df.sum(axis=1), axis=0)
    print('freq matrix\n', freq_matrix)
    ic_matrix = freq_matrix * np.log2(freq_matrix)
    print('ic matrix\n', ic_matrix)
    ic_position = ic_matrix.sum(axis=1)
    print('ic position\n', ic_position)
    # normalize the information content at each position
    norm_ic = np.log2(4) + ic_position
    print('norm ic\n', norm_ic)
    ic_matrix = freq_matrix.multiply(norm_ic, axis=0)
    print('ic matrix\n', ic_matrix)
    # calculate the information content at each position
    logo = lm.Logo(ic_matrix, color_scheme='classic',
                     width=.8, stack_order='small_on_top',
                        font_name='Arial')

    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.ax.set_title('DNA Binding Motif', fontsize=16)
    logo.ax.set_ylabel('Information Content', fontsize=14)
    logo.ax.set_xlabel('Position', fontsize=14)

    plt.show()




if __name__ == '__main__':
    # ml_sample01()
    ml_sample02()