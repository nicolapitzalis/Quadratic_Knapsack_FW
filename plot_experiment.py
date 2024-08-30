import json
import numpy as np
import matplotlib.pyplot as plt

def plot_experiments(paths):

    for path in paths:
        with open(path, 'r') as f:
            data = json.load(f)

        # y log scale
        plt.yscale('log')
        plt.plot(range(len(data['geo_mean_gaps'])), data['geo_mean_gaps'], label=data['alg'])
        
    plt.legend()
    plt.show()

if __name__ == '__main__':
    paths = ['fw_og_100_0.6_1.1_0.999.json', 'fw_step_away_100_0.6_1.1_0.999.json', 'fw_pairwise_100_0.6_1.1_0.999.json']
    plot_experiments(paths)