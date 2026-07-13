# import xarray as xr
# import matplotlib.pyplot as plt

# print([v for v in xr.open_dataset("Output.nc", group="profiles").data_vars])

# # 1. Load the dataset (pointing to your thermodynamic group)
# filename = "Output.nc"
# ds = xr.open_dataset(filename, group="profiles")

# '''
# 'pressure', 'q_rim', 'q_tot', 'ρq_ice', 'q_ice', 'q_vap', 'N_liq', 'theta_dry', 
# 'N_rai', 'q_liqonice', 'N_ice', 'q_rai', 'ρq_vap', 'ρq_liqonice', 'q_liq', 'B_rim', 
# 'temperature', 'density', 'q_sno', 'N_aer', 'theta_liq_ice', 'ρq_rai', 'ρq_rim', 
# 'SN_liq_prc', 'SN_aer_prc', 'Sq_liq_prc', 'SN_rai_prc', 'Sq_rai_prc', 'SN_liq_act', 
# 'SN_aer_act']
# '''

# temperature = ds["temperature"]
# q_tot = ds["q_tot"]
# q_ice = ds["q_ice"]
# q_vap = ds["q_vap"]



# ds["temperature"].transpose("zc", "t").plot(cmap="magma", robust=True, cbar_kwargs={"label": "Temperature (K)"}
# )
# plt.title("Temperature Profile")
# plt.xlabel("Time (s)")
# plt.ylabel("Temp (K)")
# plt.show()

# # air density
# ds["density"].transpose("zc", "t").plot(cmap="viridis", robust=True, cbar_kwargs={"label": "Density (kg/m³)"}
# )
# plt.title("Air Density Profile")
# plt.ylabel("Air density")
# plt.xlabel("Time (s)")

# plt.tight_layout()
# plt.show()

# ds["q_vap"].transpose("zc", "t").plot(cmap="YlGnBu", robust=True, cbar_kwargs={"label": "q_vap (kg/kg)"}
# )
# plt.title("Water Vapor Mixing ratio")
# plt.xlabel("Time (s)")
# plt.show()


# ds["q_sno"].transpose("zc", "t").plot(cmap="YlGnBu", robust=True, cbar_kwargs={"label": "q_sno (kg/kg)"}
# )
# plt.title("Q_sno")
# plt.xlabel("Time (s)")
# plt.show()

# # air density
# ds["density"].transpose("zc", "t").plot(cmap="viridis", robust=True, cbar_kwargs={"label": "Density (kg/m³)"}
# )
# plt.title("Air Density Profile")
# plt.ylabel("Air densit")
# plt.xlabel("Time (s)")

# plt.tight_layout()
# plt.show()

# ds["q_tot"].transpose("zc", "t").plot(cmap="YlGnBu", robust=True, cbar_kwargs={"label": "q_tot (kg/kg)"}
# )
# plt.title("Q_tot")
# plt.xlabel("Time (s)")
# plt.show()

# ds["q_ice"].transpose("zc", "t").plot(cmap="YlGnBu", robust=True, cbar_kwargs={"label": "q_ice (kg/kg)"}
# )
# plt.title("ice mixing ratio")
# plt.xlabel("Time (s)")
# plt.show()

# ds["q_rai"].transpose("zc", "t").plot(cmap="YlGnBu", robust=True, cbar_kwargs={"label": "q_rai"}
# )
# plt.title("q_rain mixing ratio")
# plt.xlabel("Time (s)")
# plt.show()

# ds["q_liq"].transpose("zc", "t").plot(cmap="YlGnBu", robust=True, cbar_kwargs={"label": "q_liq"}
# )
# plt.title("q_liq")
# plt.xlabel("Time (s)")
# plt.show()

# ds["w"].transpose("zc", "t").plot(cmap="YlGnBu", robust=True, cbar_kwargs={"label": "w"}
# )
# plt.title("w")
# plt.xlabel("Time (s)")
# plt.show()

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

filename = "Output.nc"
ds = xr.open_dataset(filename, group="profiles")

ds["true_w"] = ds["w"] / ds["density"]

vars_to_plot = [
    ("temperature", "magma", "Temperature (K)"),
    ("density", "viridis", "Density (kg/m³)"),
    ("q_vap", "YlGnBu", "q_vap (kg/kg)"),
    ("q_tot", "YlGnBu", "q_tot (kg/kg)"),
    ("q_liq", "YlGnBu", "q_liq (kg/kg)"),
    ("q_ice", "YlGnBu", "q_ice (kg/kg)"),
    ("q_rai", "YlGnBu", "q_rai (kg/kg)"),
    ("q_sno", "YlGnBu", "q_sno (kg/kg)"),
    ("ρw", "RdBu_r", "ρw (kg/m^2s)"),
    ("w", "RdBu_r", "w (m/s)"),
]

ncols = 5
nrows = 2
fig, axes = plt.subplots(nrows, ncols, figsize=(24, 10))
axes = axes.flatten()

for i, (var, cmap, label) in enumerate(vars_to_plot):
    ax = axes[i]
    data = ds[var].transpose("zc", "t")
    
    if var in ("w", "ρw"):
        vmax = float(np.abs(data).quantile(0.98))
        data.plot(ax=ax, cmap=cmap, vmin=-vmax, vmax=vmax,
                  cbar_kwargs={"label": label}, add_colorbar=True)
    else:
        data.plot(ax=ax, cmap=cmap, robust=True,
                  cbar_kwargs={"label": label}, add_colorbar=True)
    
    ax.set_title(var)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("zc (m)")

plt.suptitle("KiD Simulation Output", fontsize=16, y=1.02)
plt.tight_layout()
plt.show()



