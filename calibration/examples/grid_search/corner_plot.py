import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import argparse
import json
import matplotlib.ticker as mticker

###

# This code, which is based on a similar code in CalibrateEDMF.jl, is designed to receive loss.nc 
# from grid_search.jl function and plot a corresponding "corner plot" where the diagnoal is the 
# marginal map per parameter averaged over all other parameters and the lower half of the matrix 
# shows contour plots for all 2-parameter combinations in loss.nc
# to run from command line:
# python corner_plot.py PATH_TO_LOSS (--log_scaled_data)

###

def main():
    parser = argparse.ArgumentParser(prog='corner_plot')
    parser.add_argument("file_path", help = "path to the loss file")
    parser.add_argument("-log", "--log_scaled_data", help = "show data in log scale", action = "store_true")
    parser.add_argument("-ss", "--shared_scale", help = "show contours with scale shared", action = "store_true")
    args = parser.parse_args()
    file_path = args.file_path
    log_scaled_data = args.log_scaled_data
    shared_scale = args.shared_scale
    symbols_file = open('symbols.txt').read()
    name_dict = json.loads(symbols_file)

    data = nc.Dataset(file_path, 'r')
    param1 = []
    param2 = []
    group_names = []
    for group in data.groups:
        param1_, param2_ = group.split('.')
        param1.append(param1_)
        param2.append(param2_)
        group_names.append(group)

    param_names = list(set(param1 + param2))

    M = np.size(param_names)
    inx_mat = np.reshape(np.linspace(1,M*M,M*M),(M,M)).astype(int)
    xl, yl = np.tril_indices_from(inx_mat, k=-1)

    group_names_rev = np.flipud(group_names)
    Z_matrix = ["" for x in range(M*M)]
    X_matrix = ["" for x in range(M*M)]
    Y_matrix = ["" for x in range(M*M)]
    for i in range(0,np.size(xl)):
        Z_matrix[inx_mat[xl[i],yl[i]]-1] = group_names_rev[i]
        y_, x_ = group_names_rev[i].split('.')
        X_matrix[inx_mat[xl[i],yl[i]]-1] = x_
        Y_matrix[inx_mat[xl[i],yl[i]]-1] = y_

    x_matrix = np.reshape(X_matrix, (M,M))
    y_matrix = np.reshape(Y_matrix, (M,M))
    z_matrix = np.reshape(Z_matrix, (M,M))
    y_, x_ = z_matrix[-1,0].split('.')
    z_matrix[M-1,M-1] = y_
    for k in range(0,M-1):
        y_, x_ = z_matrix[-1,k].split('.')
        z_matrix[k,k] = x_
        x_matrix[k,k] = z_matrix[k,k]
    x_matrix[-1,-1] = z_matrix[-1,-1]

    fig, axes = plt.subplots(M, M)
    max_z = 0.0
    min_z = 1e16
    for i in range(0,M):
        for j in range(0,M):
            labelx = name_dict.get(x_matrix[i,j])
            labely = name_dict.get(y_matrix[i,j])
            if i==j:
                # if on the diagonal: collect all fields that has this parameter,
                # take thier average with respect to the other parameter and plot
                # the marginal loss curve on log-scale by scanning vertically
                # and horiozntally from this diagonal spot

                for k in range(0,j): # scan horizontally
                    group_name = z_matrix[i,k]
                    if k==0:
                        x = np.array(data.groups[group_name].variables[x_matrix[i,j]])
                        z_diag = np.zeros_like(x)
                    z = np.squeeze(np.array(data.groups[group_name].variables["loss_data"])[:,:])
                    if log_scaled_data: z = np.log10(z)
                    max_z = np.max((max_z,np.nanmax(z)))
                    min_z = np.min((min_z,np.nanmin(z)))
                    z_diag = np.add(z_diag,np.nanmean(z, axis = 0))
                for k in range(i+1,M): # scan vertically
                    group_name = z_matrix[k,j]
                    if j==0:
                        x = np.array(data.groups[group_name].variables[x_matrix[i,j]])
                        z_diag = np.zeros_like(x)
                    z = np.squeeze(np.array(data.groups[group_name].variables["loss_data"])[:,:])
                    if log_scaled_data: z = np.log10(z)
                    max_z = np.max((max_z,np.nanmax(z)))
                    min_z = np.min((min_z,np.nanmin(z)))
                    z_diag = np.add(z_diag, np.nanmean(z, axis = 1))
                ax = axes[i][j]
                pt = ax.plot(x, z_diag)
                if i==M-1:
                    ax.set_xlabel(labelx)
                else:
                    xlabels = [item.get_text() for item in ax.get_xticklabels()]
                    xempty_string_labels = [''] * len(xlabels)
                    ax.xaxis.set_major_locator(mticker.FixedLocator(ax.get_xticks().tolist()))
                    ax.set_xticklabels(xempty_string_labels)
                ax.yaxis.tick_right()

            # if off diagonal, identify and load relevant 2D slice for contour map
            elif bool(z_matrix[i,j]):
                group_name = z_matrix[i,j]
                x = np.array(data.groups[group_name].variables[x_matrix[i,j]])
                y = np.array(data.groups[group_name].variables[y_matrix[i,j]])
                z = np.squeeze(np.array(data.groups[group_name].variables["loss_data"])[:,:])
                if log_scaled_data: z = np.log10(z)
                ax = axes[i][j]
                _levels = np.linspace(min_z, max_z, 50) if shared_scale else 50
                pcm = ax.contourf(x, y, np.fliplr(np.rot90(z, k=3)), cmap = "RdYlBu_r", levels=_levels)
                if j==0 and i==M-1:
                    ax.set_xlabel(labelx)
                    ax.set_ylabel(labely)
                elif j==0:
                    xlabels = [item.get_text() for item in ax.get_xticklabels()]
                    xempty_string_labels = [''] * len(xlabels)
                    ax.xaxis.set_major_locator(mticker.FixedLocator(ax.get_xticks().tolist()))
                    ax.set_xticklabels(xempty_string_labels)
                    ax.set_ylabel(labely)
                elif i==M-1:
                    ylabels = [item.get_text() for item in ax.get_yticklabels()]
                    yempty_string_labels = [''] * len(ylabels)
                    ax.yaxis.set_major_locator(mticker.FixedLocator(ax.get_yticks().tolist()))
                    ax.set_yticklabels(yempty_string_labels)
                    ax.set_xlabel(labelx)
                else:
                    xlabels = [item.get_text() for item in ax.get_xticklabels()]
                    xempty_string_labels = [''] * len(xlabels)
                    ax.xaxis.set_major_locator(mticker.FixedLocator(ax.get_xticks().tolist()))
                    ax.set_xticklabels(xempty_string_labels)
                    ylabels = [item.get_text() for item in ax.get_yticklabels()]
                    yempty_string_labels = [''] * len(ylabels)
                    ax.yaxis.set_major_locator(mticker.FixedLocator(ax.get_yticks().tolist()))
                    ax.set_yticklabels(yempty_string_labels)

            else:
                axes[i,j].set_visible(False)
            axes[i,j].tick_params(axis='both', labelsize=8)
    fig.colorbar(pcm, ax=axes, panchor = "NE", location='top', aspect=50, anchor = (1.0,0.3))
    plt.show()

def get_loss(z_full):
    z_ens = np.mean(z_full, axis = 0)
    z_case = np.squeeze(z_ens[0,:,:])
    return z_case

if __name__ == '__main__':
    main()

