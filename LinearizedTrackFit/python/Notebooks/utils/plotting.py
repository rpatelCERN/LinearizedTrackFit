__author__ = 'demattia'

from array import array
from ROOT import TFile, TCanvas, TGraphErrors
import numpy as np
import matplotlib.pyplot as plt
from combinations import combination_index


def compute_bins(pts, pt_min, pt_max, total_bins):
    bin_size = (pt_max - pt_min)/total_bins
    pt_bins = []
    for pt in pts:
        pt_bins.append(int((pt-pt_min)/bin_size)+1)
    return pt_bins


def fit_slices(input_file_name, h_name, bin_values, value_min, value_max, total_bins,
               mean_y_min, mean_y_max, sigma_y_min, sigma_y_max,
               x_title, y_title, sigma_y_title):
    input_file = TFile(input_file_name, "READ")
    h = input_file.FindObjectAny(h_name)
    h.SetTitle("")
    c1 = TCanvas("c1", "c1", 600, 600)
    c1.SetRightMargin(0.15)
    c1.SetBottomMargin(0.15)
    h.Draw()
    c1.Update()
    h.Draw("COLZ")
    c1.Update()
    h.Draw("COLZ")
    h.GetXaxis().SetRangeUser(bin_values[0], bin_values[-1])
    h.SetTitle("")
    h.GetXaxis().SetTitle(x_title)
    h.GetXaxis().SetTitleOffset(1.4)
    h.GetXaxis().SetLabelSize(0.035)
    h.GetXaxis().SetTitleSize(0.035)
    h.GetYaxis().SetTitle(y_title)
    h.GetYaxis().SetTitleOffset(1.6)
    h.GetYaxis().SetTitleSize(0.035)
    h.GetYaxis().SetLabelSize(0.035)
    pal = h.GetListOfFunctions().FindObject("palette")
    if pal:
        pal.GetAxis().SetLabelSize(0.035)
        pal.GetAxis().SetTitleSize(0.035)
        pal.GetAxis().SetTitleOffset(1.4)
        pal.GetAxis().SetTitle("Entries")
        c1.Update()
    # c1

    bin_edges = compute_bins(bin_values, value_min, value_max, total_bins)
    print "bin_edges =", bin_edges
    projections = []
    bin_center = []
    bin_size = []
    mean = []
    mean_error = []
    sigma = []
    sigma_error = []

    c2 = TCanvas("c2", "c2", 400, 800)
    c2.Divide(2, len(bin_edges)/2)
    for i in range(len(bin_edges)-1):
        c2.cd(i+1)
        bin_center.append((bin_values[i+1]+bin_values[i])/2.)
        bin_size.append((bin_values[i+1]-bin_values[i])/2.)
        proj_y = h.ProjectionY(h.GetName()+"_py"+str(i), bin_edges[i], bin_edges[i+1]-1)
        mean.append(0.)
        mean_error.append(0.)
        sigma.append(0.)
        sigma_error.append(0.)
        if proj_y.GetEntries() > 0:
            projections.append(proj_y)
            proj_y.Fit("gaus")
            try:
                mean[-1] = proj_y.GetFunction("gaus").GetParameter(1)
                mean_error[-1] = proj_y.GetFunction("gaus").GetParError(1)
                sigma[-1] = proj_y.GetFunction("gaus").GetParameter(2)
                sigma_error[-1] = proj_y.GetFunction("gaus").GetParError(2)
            except:
                pass

    a_bin_center = array('d', bin_center)
    a_bin_size = array('d', bin_size)

    h_means = TGraphErrors(len(bin_edges), a_bin_center, array('d', mean), a_bin_size, array('d', mean_error))
    h_sigmas = TGraphErrors(len(bin_edges), a_bin_center, array('d', sigma), a_bin_size, array('d', sigma_error))

    c3 = TCanvas("Resolutions", "Resolutions", 800, 400)
    c3.Divide(2, 1)
    c3.cd(1)
    h_means.Draw("APE")
    h_means.SetTitle("")
    h_means.SetMarkerStyle(20)
    h_means.SetMarkerSize(0.5)
    h_means.GetYaxis().SetRangeUser(mean_y_min, mean_y_max)
    h_means.GetXaxis().SetTitle(x_title)
    h_means.GetYaxis().SetTitleOffset(1.5)
    h_means.GetYaxis().SetTitle(y_title + " means")
    c3.cd(2)
    h_sigmas.Draw("APE")
    h_sigmas.SetTitle("")
    h_sigmas.SetMarkerStyle(20)
    h_sigmas.SetMarkerSize(0.5)
    h_sigmas.GetYaxis().SetRangeUser(sigma_y_min, sigma_y_max)
    h_sigmas.GetXaxis().SetTitle(x_title)
    h_sigmas.GetYaxis().SetTitleOffset(1.5)
    h_sigmas.GetYaxis().SetTitle(sigma_y_title)

    return c1, c2, c3, h_means, h_sigmas, projections, h, bin_values, value_min, value_max, total_bins, input_file


def draw_chi2ndof(input_file_name):
    input_file = TFile(input_file_name, "READ")
    h = input_file.FindObjectAny("normChi2_summary")
    h.SetTitle("")
    print h
    c = TCanvas("c_chi2Ndof", "c_chi2Ndof", 600, 600)
    h.Draw()
    h.GetXaxis().SetTitle("\chi^2/ndof")
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetTitle("Entries")
    h.GetYaxis().SetTitleOffset(1.4)
    return c, h, input_file


def gaussian_fit(input_file_name, h_name):
    """
    Perform a gaussian fit on the histogram with the given h_name in the file of the given input_file_name.
    Returns the histogram, an array with the fit results (mean, mean error, sigma, sigma error) and the file.
    """
    input_file = TFile(input_file_name, "READ")
    h = input_file.FindObjectAny(h_name)
    h.Fit("gaus")
    fit_results = [0]*4
    try:
        fit_results[0] = h.GetFunction("gaus").GetParameter(1)
        fit_results[1] = h.GetFunction("gaus").GetParError(1)
        fit_results[2] = h.GetFunction("gaus").GetParameter(2)
        fit_results[3] = h.GetFunction("gaus").GetParError(2)
    except:
        pass
    h.SetName(input_file_name)
    return h, fit_results, input_file


def fit_and_draw(base_var, h_name, layers_list, radius_list):
    comb_index = combination_index(layers_list, radius_list)
    h = [gaussian_fit(base_var+"_"+str(comb_index)+".root", h_name.replace("summary", str(comb_index)))]
    for layer in layers_list:
        layers_list_removed = list(layers_list)
        radius_list_removed = list(radius_list)
        layers_list_removed.remove(layer)
        del radius_list_removed[layers_list.index(layer)]
        comb_index = combination_index(layers_list_removed, radius_list_removed)
        h.append(gaussian_fit(base_var+"_"+str(comb_index)+".root", h_name.replace("summary", str(comb_index))))
    fit_results = []
    for i in range(len(h)):
        fit_results.append(h[i][1])
    return h, fit_results


def fit_all(base_vars, h_name, layers_list, radius_list):
    h = []
    fit_results = []
    for base_var in base_vars:
        h_temp, f_temp = fit_and_draw(base_var, h_name, layers_list, radius_list)
        h.append(h_temp)
        fit_results.append(f_temp)
    return h, fit_results


def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()


def plot_combinations(fit_results, index, y_label, layers_list, legend_names, y_min=0, y_max=0):
    N = 7
    value = []
    std = []
    for f in fit_results:
        value.append([v[index] for v in f])
        std.append([v[index+1] for v in f])
    bottom_means = [0]*N
    width = 1.2       # the width of the bars
    colors = ['r', 'g', 'b', 'y']
    ind = []
    rects = []
    fig, ax = plt.subplots()
    for i in range(len(fit_results)):
        ind.append(np.linspace(i*width, 6*N+i*width, N))  # the x locations for the groups
        rects.append(ax.bar(ind[i], value[i], width, bottom=bottom_means, color=colors[i], yerr=std[i]))

    # add some text for labels, title and axes ticks
    ax.set_ylabel(y_label)
    # ax.set_title('pT relative bias')
    ax.set_xticks(ind[int(len(value)/2)])
    x_tick_labels = ['All']
    for layer in layers_list:
        x_tick_labels.append('No layer '+str(layer))
    ax.set_xticklabels(x_tick_labels)
    ax.grid()
    if y_min != y_max:
        ax.set_ylim([y_min, y_max])

    autolabel(rects[int(len(rects)/2)])

    ax.legend(rects, legend_names)

    return plt