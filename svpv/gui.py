# # -*- coding: utf-8 -*-
# """
# author: Jacob Munro, Victor Chang Cardiac Research Institute
# """
from __future__ import print_function
try:
    import Tkinter as tk
    import tkFileDialog
except ImportError:
    import tkinter as tk
    from tkinter import filedialog as tkFileDialog
from . import gui_widgets as gw
from .plot import Plot


class SVPVGui(tk.Tk):
    def __init__(self, par):
        tk.Tk.__init__(self)
        self.par = par
        self.option_add("*Font", "arial 10")
        self.buttons_1 = None
        self.reset = None
        self.list = None
        self.buttons_2 = None
        self.viewGTs = None
        self.viewSV = None
        self.plotAll = None
        self.display_cb = None
        self.setup_static_features()
        self.config(menu=gw.MenuBar(self))
        self.current_samples = []
        self.svs = self.par.run.vcf.filter_svs(self.par.filter)
        self.sample_selector = None
        self.set_sample_selector()
        self.genotype_selector = None
        self.set_genotype_selector()
        self.plot_custom = None
        self.set_plot_custom()
        self.filters = None
        self.set_filters()
        self.sv_chooser = None
        self.set_sv_chooser()
        self.info_box = None
        self.set_info_box()
        self.filename = None

    def setup_static_features(self):
        self.wm_title("SVPV - Structural Variant Prediction Viewer")
        self.window_size()

        if self.buttons_1:
            self.buttons_1.destroy()
        self.buttons_1 = tk.LabelFrame(self)
        if self.reset:
            self.reset.destroy()
        self.reset = tk.Button(self.buttons_1, text="Reset Filters", command=self.reset_filters)
        self.reset.grid(row=0, column=0, padx=40, sticky=tk.W)
        if self.list:
            self.list.destroy()
        self.list = tk.Button(self.buttons_1, text="Apply Filters", command=self.apply_filters)
        self.list.grid(row=0, column=1, padx=40, sticky=tk.E)
        self.buttons_1.grid(row=4, column=0, columnspan=2, sticky=tk.EW, padx=10)

        if self.buttons_2:
            self.buttons_2.destroy()
        self.buttons_2 = tk.LabelFrame(self)
        if self.viewGTs:
            self.viewGTs.destroy()
        self.viewGTs = tk.Button(self.buttons_2, text="Get Genotypes", command=self.view_gts)
        self.viewGTs.grid(row=0, column=0, padx=25, sticky=tk.W)
        if self.viewSV:
            self.viewSV.destroy()
        self.viewSV = tk.Button(self.buttons_2, text="Plot Selected SV", command=self.plot_sv)
        self.viewSV.grid(row=0, column=1, padx=25, sticky=tk.W)
        if self.plotAll:
            self.plotAll.destroy()
        self.plotAll = tk.Button(self.buttons_2, text="Plot All SVs", command=self.plot_all)
        self.plotAll.grid(row=0, column=2, padx=25, sticky=tk.E)
        if self.display_cb:
            self.display_cb.destroy()
        self.display_var = tk.IntVar(value=1)
        self.display_cb = tk.Checkbutton(self.buttons_2, text='display plot on creation', variable=self.display_var,
                                         onvalue=1, offvalue=0)
        self.display_cb.grid(row=1, column=1, padx=25, sticky=tk.E)
        self.buttons_2.grid(row=6, column=0, columnspan=2, sticky=tk.EW, padx=10)

    def text_size(self, opt):
        if opt == 1:
            self.option_add("*Font", "arial 8")
        elif opt == 2:
            self.option_add("*Font", "arial 10")
        elif opt == 3:
            self.option_add("*Font", "arial 12")
        elif opt == 4:
            self.option_add("*Font", "arial 14")
        elif opt == 5:
            self.option_add("*Font", "arial 16")
        self.setup_static_features()
        self.set_sample_selector()
        self.set_genotype_selector()
        self.set_filters()
        self.set_sv_chooser()
        self.set_info_box()

    def set_sample_selector(self):
        if self.sample_selector:
            self.sample_selector.destroy()
        self.sample_selector = gw.SampleSelector(self, self.par.run.samples)
        self.sample_selector.grid(row=1, column=0, sticky=tk.NSEW, padx=10)

    def set_genotype_selector(self):
        if self.genotype_selector:
            self.genotype_selector.destroy()
        self.genotype_selector = gw.SampleGenotypeSelector(self, self.current_samples)
        self.genotype_selector.grid(row=1, column=1, sticky=tk.NSEW, padx=10)

    def set_plot_custom(self):
        if self.plot_custom:
            self.plot_custom.destroy()
        self.plot_custom = gw.PlotCustom(self)
        self.plot_custom.grid(row=2, column=0, columnspan=2, sticky=tk.NSEW, pady=2, padx=10)

    def set_filters(self):
        if self.filters:
            self.filters.destroy()
        self.filters = gw.Filters(self)
        self.filters.grid(row=3, column=0, columnspan=2, sticky=tk.NSEW, pady=2, padx=10)

    def set_sv_chooser(self):
        if self.sv_chooser:
            self.sv_chooser.destroy()
        self.sv_chooser = gw.SvChooser(self, self.svs, self.par.run.vcf.count)
        self.sv_chooser.grid(row=5, column=0, sticky=tk.NSEW, padx=10, columnspan=2)

    def set_info_box(self, message=''):
        if self.info_box:
            self.info_box.destroy()
        self.info_box = gw.InfoBox(self, message)
        self.info_box.grid(row=7, column=0, sticky=tk.NSEW, padx=10, columnspan=2)

    def reset_filters(self):
        self.set_info_box()
        self.current_samples = []
        self.set_genotype_selector()
        self.filters.reset()
        self.apply_filters()
        self.set_sv_chooser()

    def switch_vcf(self, i):
        print('switching to {}'.format(self.par.run.alt_vcfs[i].name))
        tmp = self.par.run.vcf
        self.par.run.vcf = self.par.run.alt_vcfs[i]
        self.par.run.alt_vcfs[i] = tmp
        self.config(menu=gw.MenuBar(self))
        self.reset_filters()

    def view_gts(self):
        self.set_info_box()
        self.info_box.genotypes(self.svs[self.sv_chooser.sv_fl.selected_idx], self.par.run.vcf.samples, self.par.run.samples)

    def apply_filters(self):
        self.set_info_box()
        self.par.filter.sample_GTs = {}
        if self.genotype_selector.GT_CBs:
            for i, gt_cb in enumerate(self.genotype_selector.GT_CBs):
                if gt_cb.get_selection():
                    self.par.filter.sample_GTs[self.genotype_selector.samples[i]] = gt_cb.get_selection()
        else:
            self.par.filter.GTs = ['*']

        # update af filter
        self.par.filter.AF_thresh = None
        if self.filters.af_filter.af_on.get():
            self.par.filter.AF_thresh = self.filters.af_filter.af_var.get()
            if self.filters.af_filter.af_gt_lt.get() == 1:
                self.par.filter.AF_thresh_is_LT = True
            else:
                self.par.filter.AF_thresh_is_LT = False

        # update gene_list_intersection
        self.par.filter.gene_list_intersection = False
        if self.filters.gene_filter.gene_list_on.get():
            self.par.filter.gene_list_intersection = True

        # update ref_gene_intersection
        self.par.filter.RG_intersection = False
        if self.filters.gene_filter.ref_gene_on.get():
            self.par.filter.RG_intersection = True

        # update exonic
        self.par.filter.exonic = False
        if self.filters.gene_filter.exonic_on.get():
            self.par.filter.exonic = True

        # update length filter
        self.par.filter.max_len = None
        self.par.filter.min_len = None
        if self.filters.len_filter.len_GT_On.get():
            units = 1
            if str(self.filters.len_filter.len_GT_Units.get()) == 'kbp':
                units = 1000
            elif str(self.filters.len_filter.len_GT_Units.get()) == 'Mbp':
                units = 1000000
            self.par.filter.min_len = units * int(self.filters.len_filter.len_GT_val.get())
        if self.filters.len_filter.len_LT_On.get():
            units = 1
            if str(self.filters.len_filter.len_LT_Units.get()) == 'kbp':
                units = 1000
            elif str(self.filters.len_filter.len_LT_Units.get()) == 'Mbp':
                units = 1000000
            self.par.filter.max_len = units * int(self.filters.len_filter.len_LT_val.get())

        # update svtype filter
        self.par.filter.svtype = gw.SvTypeFilter.types[self.filters.type_filter.type_var.get()]

        # get list of svs based on current filters
        self.svs = self.par.run.vcf.filter_svs(self.par.filter)
        self.set_sv_chooser()

    def plot_sv(self, sv=None):
        self.set_info_box()
        if not self.current_samples:
            self.info_box.message.config(text="Error: No Samples Selected")
        elif self.sv_chooser.sv_fl.selected_idx is None and sv is None:
            self.info_box.message.config(text="Error: No SV Selected")
        else:
            if not sv:
                sv = self.svs[self.sv_chooser.sv_fl.selected_idx]
            plot = Plot(sv, self.current_samples, self.par)
            if self.display_var.get():
                self.filename = plot.plot_figure(group=self.par.plot.grouping, display=self.par.run.display)
            else:
                self.filename = plot.plot_figure(group=self.par.plot.grouping, display=False)
            self.info_box.message.config(text='plot path copied to clipboard')
            self.clipboard_clear()
            self.clipboard_append(self.filename)
            self.update()

    def set_plot_all_dir(self):
        dir_options = {}
        dir_options['initialdir'] = self.par.run.out_dir
        dir_options['parent'] = self
        dir_options['title'] = 'select existing or type new directory'
        dir_options['mustexist'] = False
        path = tkFileDialog.askdirectory(**dir_options)
        if path == '':
            return None
        else:
            return path

    def plot_all(self):
        self.set_info_box()
        if not self.current_samples:
            self.info_box.message.config(text="Error: No Samples Selected")
        else:
            old_path = self.par.run.out_dir
            new_path = self.set_plot_all_dir()
            if not new_path:
                return None
            self.par.run.out_dir = new_path
            for i, sv in enumerate(self.svs):
                message = "Plotting %d of %d." % (i+1, len(self.svs))
                self.info_box.message.config(text=message)
                self.update()
                plot = Plot(sv, self.current_samples, self.par)
                self.filename = plot.plot_figure()
            self.info_box.message.config(text="Done.")
            self.update()
            self.par.run.out_dir = old_path

    def window_size(self):
        sw = self.winfo_screenwidth()
        sh = self.winfo_screenheight()
        x = int(0.05 * sw)
        y = int(0.05 * sh)
        self.geometry('+%d+%d' % (x, y))

    def samples_update(self, idxs):
        self.set_info_box()
        self.current_samples = []
        for idx in idxs:
            self.current_samples.append(self.par.run.samples[int(idx)])
        self.genotype_selector.destroy()
        self.genotype_selector = gw.SampleGenotypeSelector(self, self.current_samples)
        self.genotype_selector.grid(row=1, column=1, sticky=tk.NSEW, padx=10)


def main(par):
    root = SVPVGui(par)
    root.mainloop()

