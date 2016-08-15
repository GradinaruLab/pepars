
def plot_amino_acid_property_distribution_from_matrix(matrix_of_interest,property_name_string,plot_title):
    sequence_length = matrix_of_interest.shape[1]
    max_of_matrix = matrix_of_interest.max()
    min_of_matrix = matrix_of_interest.min()
    plt.title(plot_title)
    for i,v in enumerate(xrange(sequence_length)):
        v = v + 1
        ax1 = subplot(sequence_length,1,v)
        ax1.set_xlim(min_of_matrix,max_of_matrix)
        sns.distplot(matrix_of_interest[:,i],bins=100,kde=False,rug=False,ax = ax1)
    plt.subplots_adjust(hspace=0.5)
    plt.xlabel(property_name_string)
    plt.show()



def plot_amino_acid_property_distribution_from_array(array_of_interest_high,array_of_interest_low,property_name_string,plot_title_high,plot_title_low):
    # Determine minimum and maximum for the axis 
    max_of_array_high = max(array_of_interest_high)
    max_of_array_low =  max(array_of_interest_low)
    min_of_array_high = min(array_of_interest_high)
    min_of_array_low = min(array_of_interest_low)
    max_array = max(max_of_array_high,max_of_array_low)
    min_array = min(min_of_array_high,min_of_array_low)
    f,(ax1,ax2) = plt.subplots(2)
    sns.distplot(array_of_interest_high,bins=100,kde=False,ax=ax1)
    plt.xlabel(property_name_string)
    plt.ylabel('sequence counts')
    plt.title(plot_title_high)
    sns.distplot(array_of_interest_low,bins=100,kde=False,ax=ax2)
    plt.xlabel(property_name_string)
    plt.ylabel('sequence counts')
    plt.title(plot_title_low)
    plt.show()