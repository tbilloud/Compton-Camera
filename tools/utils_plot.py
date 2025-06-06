def plot_hitsNclusters(pixelHits, pixelClusters, max_keV):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(7, 8))

    axes[0].hist(pixelHits['Energy (keV)'], bins=max_keV, range=(0, max_keV))
    axes[0].set_xlabel('Energy (keV)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Pixel Hits')

    axes[1].hist(pixelClusters['Energy (keV)'], bins=max_keV, range=(0, max_keV))
    axes[1].set_xlabel('Energy (keV)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Single Clusters')

    plt.tight_layout()
    plt.show()