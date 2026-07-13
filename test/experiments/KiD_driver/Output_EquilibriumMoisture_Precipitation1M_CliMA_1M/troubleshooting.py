import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

filename = "Output.nc"
ds = xr.open_dataset(filename, group="profiles")

vars_to_plot = [
    ("temperature", "magma", "Temperature (K)"),
    ("density", "viridis", "Density (kg/m³)"),
    ("q_vap", "YlGnBu", "q_vap (kg/kg)"),
    ("q_tot", "YlGnBu", "q_tot (kg/kg)"),
    ("q_liq", "YlGnBu", "q_liq (kg/kg)"),
    ("q_ice", "YlGnBu", "q_ice (kg/kg)"),
    ("q_rai", "YlGnBu", "q_rai (kg/kg)"),
    ("q_sno", "YlGnBu", "q_sno (kg/kg)"),
    ("w", "RdBu_r", "w (m/s)"),
    ("w", "RdBu_r", "w (m/s)"),
]

ncols = 5
nrows = 2
fig, axes = plt.subplots(nrows, ncols, figsize=(24, 10))
axes = axes.flatten()

for i, (var, cmap, label) in enumerate(vars_to_plot):
    ax = axes[i]
    data = ds[var].transpose("zc", "t")
    
    # Center colormap at 0 for w
    if var == "w":
        vmax = float(np.abs(data).quantile(0.98))
        im = data.plot(ax=ax, cmap=cmap, vmin=-vmax, vmax=vmax,
                      cbar_kwargs={"label": label}, add_colorbar=True)
    else:
        im = data.plot(ax=ax, cmap=cmap, robust=True,
                      cbar_kwargs={"label": label}, add_colorbar=True)
    
    ax.set_title(var)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("zc (m)")

plt.suptitle("KiD Simulation Output", fontsize=16, y=1.02)
plt.tight_layout()
plt.show()