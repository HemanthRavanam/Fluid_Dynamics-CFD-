import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV file
df = pd.read_csv('CFD2-LAB2.1-Ravanam-Lonkar-Meesala-plot-data.csv')
x = np.linspace(0, 1, 102)

titles = [
    "Velocity U at x = 0.8 Re 1", "Velocity V at y = 0.6 Re 1", "Pressure at y = 0.6  Re 1",
    "Velocity U at x = 0.8 Re 100", "Velocity V at y = 0.6 Re 100", "Pressure at y = 0.6 Re 100",
    "Velocity U at x = 0.8 Re 1000", "Velocity V at y = 0.6 Re 1000", "Pressure at y = 0.6 Re 1000"
]

# Plotting
fig, axes = plt.subplots(3, 3, figsize=(15, 12))
axes = axes.flatten()

for i in range(len(df.columns) // 2):
    col1 = df.columns[2*i]
    col2 = df.columns[2*i + 1]
    ax = axes[i]
    
    ax.plot(x, df[col1], label=col1)
    ax.plot(x, df[col2], label=col2)
    
    ax.set_title(titles[i], fontsize=10)
    ax.set_xlabel('along line')
    ax.set_ylabel('Values')
    ax.grid(True)
    ax.legend()

plt.tight_layout()
plt.show()