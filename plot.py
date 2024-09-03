import json
import numpy as np
import matplotlib.pyplot as plt

def plot_experiments(paths):

    for path in paths:
        with open(path, 'r') as f:
            data = json.load(f)

        # y log scale
        plt.yscale('log')
        plt.plot(range(len(data['mean_gaps'])), data['mean_gaps'], label=data['alg'])
        
    plt.legend()
    plt.show()

if __name__ == '__main__':
    p = "/Users/lorenzoleuzzi/Documents/MATLAB/"
    paths = [f'{p}/fw_og_100_0.6_1_0.001.json', f'{p}/fw_away_step_100_0.6_1_0.001.json', f'{p}/fw_pairwise_100_0.6_1_0.001.json']
    # paths = [f'{p}/fw_pairwise_100_0.6_1.1_1.json']
    plot_experiments(paths)