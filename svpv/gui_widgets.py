# # -*- coding: utf-8 -*-
# """
# author: Jacob Munro, Victor Chang Cardiac Research Institute
# """
from __future__ import print_function
import Tkinter as tk
import tkFileDialog
import re
from shutil import copyfile


class MenuBar(tk.Menu):
    def __init__(self, parent):
        tk.Menu.__init__(self, parent)
        self.parent = parent

        self.file = tk.Menu(self, tearoff=0)
        self.file.add_command(label='save figure as', command=self.asksaveasfilename)
        self.file.add_command(label='exit', command=self.parent.quit)
        self.add_cascade(label='File', menu=self.file)

        self.size_var = tk.IntVar(value=2)
        self.resize = tk.Menu(self, tearoff=0)
        self.resize.add_radiobutton(label='x-small', var=self.size_var, value=1, command=self.update_size)
        self.resize.add_radiobutton(label='small', var=self.size_var, value=2, command=self.update_size)
        self.resize.add_radiobutton(label='medium', var=self.size_var, value=3, command=self.update_size)
        self.resize.add_radiobutton(label='large', var=self.size_var, value=4, command=self.update_size)
        self.resize.add_radiobutton(label='x-large', var=self.size_var, value=5, command=self.update_size)
        self.add_cascade(label='View', menu=self.resize)

        self.vcf_var = tk.IntVar(value=0)
        self.vcf = tk.Menu(self, tearoff=0)
        self.vcf.add_radiobutton(label=parent.par.run.vcf.name, var=self.vcf_var, value=0,
                                        command=self.switch_vcf)
        for i in range(len(parent.par.run.alt_vcfs)):
            self.vcf.add_radiobutton(label=parent.par.run.alt_vcfs[i].name, var=self.vcf_var, value=i+1,
                                            command=self.switch_vcf)
        self.add_cascade(label='VCF', menu=self.vcf)

    def asksaveasfilename(self):
        if not self.parent.filename:
            self.parent.set_info_box(message='Error: No figure has been created yet.')
        else:
            file_options = {}
            file_options['initialdir'] = self.parent.par.run.out_dir
            file_options['initialfile'] = re.sub('/.+/', '', self.parent.filename)
            file_options['filetypes'] = [('pdf files', '.pdf')]
            file_options['parent'] = self.parent
            file_options['title'] = 'save figure as'
            filename = tkFileDialog.askdirectory(**file_options)
            if filename:
                copyfile(self.parent.filename, filename)

    def update_size(self):
        self.parent.text_size(self.size_var.get())

    def switch_vcf(self):
        if self.vcf_var.get() != 0:
            self.parent.switch_vcf(self.vcf_var.get()-1)


class SampleSelector(tk.LabelFrame):
    def __init__(self, parent, samples):
        tk.LabelFrame.__init__(self, parent, text="Sample Selection")
        self.parent = parent
        self.samples = samples
        self.lb = tk.Listbox(self, selectmode=tk.EXTENDED, width=20)
        for s in self.samples:
            self.lb.insert(tk.END, s)
        self.lb.grid(row=0, column=0, padx=10, columnspan=2)
        self.scroll = tk.Scrollbar(self, orient=tk.VERTICAL)
        self.scroll.grid(row=0, column=2, sticky=tk.NS)
        self.lb.config(yscrollcommand=self.scroll.set)
        self.scroll.config(command=self.lb.yview)
        self.setter = tk.Button(self, text="Select", command=self.select)
        self.setter.grid(row=1,column=0, padx=10)
        self.clear = tk.Button(self, text="Clear", command=self.clear)
        self.clear.grid(row=1, column=1, padx=10)

    def select(self):
        self.parent.samples_update(self.lb.curselection())

    def clear(self):
        self.parent.current_samples = []
        self.parent.set_genotype_selector()


class SampleGenotypeSelector(tk.LabelFrame):
    def __init__(self, parent, samples):
        tk.LabelFrame.__init__(self, parent, text="Sample Genotype Selection")
        self.parent = parent
        self.samples = samples
        self.max_row = 5
        if self.samples:
            self.GT_CBs = []
            self.c = 0
            self.r = 0
            self.set_samples()
        else:
            self.GT_CBs = None
            self.lab = tk.Label(self,text="-- No samples Selected --")
            self.lab.grid(row=0, sticky = tk.EW)

    def set_samples(self):
        for s in self.samples:
            self.GT_CBs.append(GenotypeSelector(self, s))
            self.GT_CBs[-1].grid(column=self.c, row=self.r, sticky = tk.EW)
            self.r += 1
            if self.r == self.max_row:
                self.r = 0
                self.c +=1


class GenotypeSelector(tk.LabelFrame):
    def __init__(self, parent, name):
        tk.LabelFrame.__init__(self, parent, text=name)
        self.gts = ["0/0", "0/1", "1/1"]
        self.checkVars = []
        self.CBs = []
        self.r = 0
        self.c = 0
        self.set_check_boxes()

    def set_check_boxes(self):
        for gt in self.gts:
            self.checkVars.append(tk.IntVar(value=0))
            self.CBs.append(tk.Checkbutton(self, text=gt, variable=self.checkVars[-1], onvalue=1, offvalue=0))
            self.CBs[-1].grid(row = self.r, column=self.c, sticky=tk.EW)
            self.c += 1

    def get_selection(self):
        gts = []
        for i in range(0, len(self.gts)):
            if self.checkVars[i].get():
                gts.append(self.gts[i])
        return gts


class Filters(tk.LabelFrame):
    def __init__(self, parent):
        tk.LabelFrame.__init__(self, parent, text="Filters")
        self.af_filter = AfFilter(self)
        self.gene_filter = GeneFilter(self)
        self.len_filter = LengthFilter(self)
        self.type_filter = SvTypeFilter(self)

        self.af_filter.grid(row=0, sticky=tk.NSEW, columnspan=2)
        self.gene_filter.grid(row=0, column=2, sticky=tk.NSEW)
        self.len_filter.grid(row=0, column=3, sticky=tk.NSEW)
        self.type_filter.grid(row=0, column=4, sticky=tk.NSEW)

    def reset(self):
        self.af_filter.reset()
        self.gene_filter.reset()
        self.len_filter.reset()
        self.type_filter.reset()

class SvTypeFilter(tk.LabelFrame):
    types = (None, 'DEL', 'DUP', 'CNV', 'INV', 'INS', 'BND', 'TRA')
    def __init__(self, parent):
        tk.LabelFrame.__init__(self, parent)
        self.type_var = tk.IntVar(value=0)
        self.radio_buttons = [tk.Radiobutton(self, text='All', justify=tk.LEFT, variable=self.type_var, value=0)]
        self.radio_buttons[-1].grid(row=0, column=0, sticky=tk.W)
        r=1
        c=0
        for i in range(1, len(SvTypeFilter.types)):
            self.radio_buttons.append(tk.Radiobutton(self, text=SvTypeFilter.types[i], justify=tk.LEFT,
                                                     variable=self.type_var, value=i))
            self.radio_buttons[-1].grid(row=r, column=c, sticky=tk.W)
            r += 1
            if r > 3:
                r = 0
                c += 1

    def reset(self):
        self.type_var.set(0)

class GeneFilter(tk.LabelFrame):
    def __init__(self, parent):
        tk.LabelFrame.__init__(self, parent)
        self.ref_gene_on = tk.IntVar(value=0)
        self.ref_gene_on_cb = tk.Checkbutton(self, text="ref gene\nintersection", justify=tk.LEFT, variable=self.ref_gene_on)
        self.ref_gene_on_cb.grid(row=0, column=0, sticky=tk.W)

        self.gene_list_on = tk.IntVar(value=0)
        self.gene_list_on_cb = tk.Checkbutton(self, text="gene list\nintersection", justify=tk.LEFT, variable=self.gene_list_on)
        self.gene_list_on_cb.grid(row=1, column=0, sticky=tk.W)

        self.exonic_on = tk.IntVar(value=0)
        self.exonic_on_cb = tk.Checkbutton(self, text="exonic", justify=tk.LEFT, variable=self.exonic_on)
        self.exonic_on_cb.grid(row=2, column=0, sticky=tk.W)

    def reset(self):
        self.ref_gene_on.set(0)
        self.gene_list_on.set(0)
        self.exonic_on.set(0)


class LengthFilter(tk.LabelFrame):
    def __init__(self, parent):
        tk.LabelFrame.__init__(self, parent)
        self.title = tk.Label(self, text='SV Length')

        self.len_GT_On = tk.IntVar(value=0)
        self.len_GT_On_CB = tk.Checkbutton(self, text=">", justify=tk.LEFT, variable=self.len_GT_On)
        self.len_GT_val = tk.Spinbox(self, values=(1,5,10,50,100,500), width=3)
        self.len_GT_Units = tk.Spinbox(self, values=("bp", "kbp", "Mbp"), width=3)

        self.len_LT_On = tk.IntVar(value=0)
        self.len_LT_On_CB = tk.Checkbutton(self, text="<", justify=tk.LEFT, variable=self.len_LT_On)
        self.len_LT_val = tk.Spinbox(self, values=(1,5,10,50,100,500), width=3)
        self.len_LT_Units = tk.Spinbox(self, values=("bp", "kbp", "Mbp"), width=3)

        self.title.grid(row=0, column=0, sticky=tk.EW, columnspan=4)
        self.len_GT_On_CB.grid(row=1, column=0, sticky=tk.EW)
        self.len_GT_val.grid(row=2, column=0, sticky=tk.EW)
        self.len_GT_Units.grid(row=2, column=1, sticky=tk.EW)
        self.len_LT_On_CB.grid(row=1, column=2, sticky=tk.EW)
        self.len_LT_val.grid(row=2, column=2, sticky=tk.EW)
        self.len_LT_Units.grid(row=2, column=3, sticky=tk.EW)

    def reset(self):
        self.len_GT_On.set(0)
        self.len_LT_On.set(0)
        self.len_GT_val.destroy()
        self.len_GT_Units.destroy()
        self.len_LT_val.destroy()
        self.len_LT_Units.destroy()
        self.len_GT_val = tk.Spinbox(self, values=(1, 5, 10, 50, 100, 500), width=3)
        self.len_GT_Units = tk.Spinbox(self, values=("bp", "kbp", "Mbp"), width=3)
        self.len_LT_val = tk.Spinbox(self, values=(1, 5, 10, 50, 100, 500), width=3)
        self.len_LT_Units = tk.Spinbox(self, values=("bp", "kbp", "Mbp"), width=3)
        self.len_GT_val.grid(row=2, column=0, sticky=tk.EW)
        self.len_GT_Units.grid(row=2, column=1, sticky=tk.EW)
        self.len_LT_val.grid(row=2, column=2, sticky=tk.EW)
        self.len_LT_Units.grid(row=2, column=3, sticky=tk.EW)


class AfFilter(tk.LabelFrame):
    def __init__(self, parent):
        tk.LabelFrame.__init__(self, parent)
        self.af_on = tk.IntVar(value=0)
        self.af_on_cb = tk.Checkbutton(self, text="AF Threshold", variable=self.af_on)
        self.af_on_cb.grid(row=0, column=0, padx=3, columnspan=2, sticky=tk.EW)
        self.af_var = tk.DoubleVar()
        self.af_gt_lt = tk.IntVar(value=1)
        self.af_lt_radio_button = tk.Radiobutton(self, text="<", variable=self.af_gt_lt, value=1)
        self.af_lt_radio_button.grid(row=1, column=0, sticky=tk.EW)
        self.af_gt_radio_button = tk.Radiobutton(self, text=">", variable=self.af_gt_lt, value=2)
        self.af_gt_radio_button.grid(row=1, column=1, sticky=tk.EW)

        self.af_scale = tk.Scale(self, variable=self.af_var, from_=float(0), to=float(1), resolution=float(0.01),
                                 orient=tk.HORIZONTAL)
        self.af_scale.grid(row=2, column=0, padx=3, sticky=tk.EW, columnspan=2)

    def reset(self):
        self.af_on.set(0)
        self.af_var.set(0)
        self.af_gt_lt.set(1)


class SvChooser(tk.LabelFrame):
    def __init__(self, parent, svs, sv_count):
        tk.LabelFrame.__init__(self, parent, text="Structural Variant Call Selection")
        self.sv_fl = None
        if not svs:
            self.lab = tk.Label(self,text="-- No matches --")
            self.lab.grid(row=0, column=0, sticky = tk.EW)
            self.num_svs_lab = tk.Label(self, text='-- of %d SVs' % sv_count)
        else:
            self.sv_fl = FieldedListbox(self, ("SV Type", "Chr A", "Pos A", "Chr B", "Pos B", "Length (bp)", "AF"))
            for sv in svs:
                self.sv_fl.push_entry(sv.string_tuple())
            self.num_svs_lab = tk.Label(self, text='%d of %d SVs' % (len(svs), sv_count))
            self.sv_fl.grid(row=0, sticky = tk.NSEW)
        self.num_svs_lab.grid(row=1, column=0, sticky=tk.EW)


class InfoBox(tk.LabelFrame):
    def __init__(self, parent, message):
        tk.LabelFrame.__init__(self, parent, text="Info")
        self.message = tk.Label(self, text=message, bg='white', width=55, justify=tk.LEFT)
        self.message.pack()

    def genotypes(self, sv, vcf_samples, run_samples):
        message = ''
        line = ''
        for i, gt in enumerate(sv.GTs):
            if '1' in gt:
                if vcf_samples[i] in run_samples:
                    this = '%s:%s' % (vcf_samples[i], gt)
                    if len(line) == 0:
                        line = this
                    elif len(this) + len(line) >= 50:
                        message += line + ',\n'
                        line = this
                    else:
                        line += ', %s' % this
        message += line
        self.message.config(text=message, bg='white', width=55, justify=tk.LEFT)


# class for multiple listboxes linked by single scroll bar
class FieldedListbox(tk.Frame):
    def __init__(self, parent, header):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.num_f = len(header)
        self.headers = []
        self.lbs = []
        self.c = 0
        self.selected_idx = None

        self.scroll = tk.Scrollbar(self, orient=tk.VERTICAL)
        for i in range(0, self.num_f):
            self.headers.append(tk.Label(self, text=header[i]))
            self.headers[-1].grid(row=0, column=self.c, sticky=tk.EW)
            self.lbs.append(tk.Listbox(self, width=10, yscrollcommand= self.yscroll, bg='white'))
            self.lbs[-1].grid(row=1, column=self.c, sticky=tk.EW)
            self.lbs[-1].bind("<<ListboxSelect>>", self.select)
            self.c += 1

        self.scroll.config(command=self.lbs[0].yview)
        self.scroll.grid(row=1, column=self.c, sticky=tk.NS)

    def select(self, val):
        idx = int(val.widget.curselection()[0])
        for lb in self.lbs:
            if self.selected_idx is not None:
                lb.itemconfig(self.selected_idx, background='white')
            lb.activate(idx)
            lb.itemconfig(idx, background='gray70', selectbackground='gray70')
        self.selected_idx = idx

    def yscroll(self, *args):
        self.scroll.set(*args)
        for i in range(0, self.num_f):
            self.lbs[i].yview_moveto(args[0])

    def push_entry(self, entry):
        if not len(entry) == self.num_f:
            return False
        for i, f in enumerate(entry):
            self.lbs[i].insert(tk.END, f)
        return True

    # return the index of selection
    def get_selection(self):
        return int(self.lbs[0].curselection())
