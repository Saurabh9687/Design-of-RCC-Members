#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from tkinter import *
from tkinter import ttk
from tkinter.scrolledtext import *
from tkinter import messagebox
from fpdf import FPDF
import math
from scipy.interpolate import interp2d
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import fpdf
import csv
import os
import platform
import subprocess
import sys

app = Tk()
app.title("BEAM DESIGN")
style = ttk.Style(app)

# INPUT DATA
clear_span_text = IntVar();
clear_span_text.set(0)
clear_span_label = Label(app, text="Clear Span", font=('normal', 10), padx=5, pady=5);
clear_span_label.grid(row=1, column=0);
clear_span_entry = Entry(app, textvariable=clear_span_text, width=30);
clear_span_entry.grid(row=1, column=1);
clear_span_unit = Label(app, text="m", font=('normal', 10), padx=5, pady=5);
clear_span_unit.grid(row=1, column=2);


width_Beam_text = IntVar()
width_Beam_text.set(0)
width_Beam_label = Label(app, text="Width of Beam", font=('normal', 10), padx=5, pady=5)
width_Beam_label.grid(row=2, column=0)
width_Beam_entry = Entry(app, textvariable=width_Beam_text, width=30)
width_Beam_entry.grid(row=2, column=1)
width_Beam_unit = Label(app, text="mm", font=('normal', 10), padx=5, pady=5);
width_Beam_unit.grid(row=2, column=2);


width_support_text = IntVar()
width_support_text.set(0)
width_support_label = Label(app, text="Width of support", font=('normal', 10), padx=5, pady=5)
width_support_label.grid(row=3, column=0)
width_support_entry = Entry(app, textvariable=width_support_text, width=30)
width_support_entry.grid(row=3, column=1)
width_support_unit = Label(app, text="mm", font=('normal', 10), padx=5, pady=5);
width_support_unit.grid(row=3, column=2);


clear_cover_text = IntVar()
clear_cover_text.set(0)
clear_cover_label = Label(app, text="Clear cover", font=('normal', 10), padx=5, pady=5)
clear_cover_label.grid(row=4, column=0)
clear_cover_entry = Entry(app, textvariable=clear_cover_text, width=30)
clear_cover_entry.grid(row=4, column=1)
clear_cover_unit = Label(app, text="mm", font=('normal', 10), padx=5, pady=5);
clear_cover_unit.grid(row=4, column=2);


overall_depth_text = IntVar()
overall_depth_text.set(0)
overall_depth_label = Label(app, text="Overall depth", font=('normal', 10), padx=5, pady=5)
overall_depth_label.grid(row=5, column=0)
overall_depth_entry = Entry(app, textvariable=overall_depth_text, width=30)
overall_depth_entry.grid(row=5, column=1)
overall_depth_unit = Label(app, text="mm", font=('normal', 10), padx=5, pady=5);
overall_depth_unit.grid(row=5, column=2);


dia_tension_reinforcement_text = IntVar()
dia_tension_reinforcement_text.set(0)
dia_tension_reinforcement_label = Label(app, text="Dia of tension reinforcement", font=('normal', 10), padx=5, pady=5)
dia_tension_reinforcement_label.grid(row=6, column=0)
dia_tension_reinforcement_entry = Entry(app, textvariable=dia_tension_reinforcement_text, width=30)
dia_tension_reinforcement_entry.grid(row=6, column=1)
dia_tension_reinforcement_unit = Label(app, text="mm", font=('normal', 10), padx=5, pady=5);
dia_tension_reinforcement_unit.grid(row=6, column=2);


dia_comp_reinforcement_text = IntVar()
dia_comp_reinforcement_text.set(0)
dia_comp_reinforcement_label = Label(app, text="Dia of compression reinforcement", font=('normal', 10), padx=5, pady=5)
dia_comp_reinforcement_label.grid(row=7, column=0)
dia_comp_reinforcement_entry = Entry(app, textvariable=dia_comp_reinforcement_text, width=30)
dia_comp_reinforcement_entry.grid(row=7, column=1)
dia_comp_reinforcement_unit = Label(app, text="mm", font=('normal', 10), padx=5, pady=5);
dia_comp_reinforcement_unit.grid(row=7, column=2);


dia_stirrups_text = IntVar()
dia_stirrups_text.set(0)
dia_stirrups_label = Label(app, text="Dia of stirrups", font=('normal', 10), padx=5, pady=5)
dia_stirrups_label.grid(row=8, column=0)
dia_stirrups_entry = Entry(app, textvariable=dia_stirrups_text, width=30)
dia_stirrups_entry.grid(row=8, column=1)
dia_stirrups_unit = Label(app, text="mm", font=('normal', 10), padx=5, pady=5);
dia_stirrups_unit.grid(row=8, column=2);


grade_stirrups_text = IntVar()
grade_stirrups_text.set(0)
grade_stirrups_label = Label(app, text="Grade of stirrups", font=('normal', 10), padx=5, pady=5)
grade_stirrups_label.grid(row=9, column=0)
grade_stirrups_entry = Entry(app, textvariable=grade_stirrups_text, width=30)
grade_stirrups_entry.grid(row=9, column=1)
grade_stirrups_unit = Label(app, text="N/mm2", font=('normal', 10), padx=5, pady=5);
grade_stirrups_unit.grid(row=9, column=2);


grade_concrete_text = IntVar()
grade_concrete_text.set(0)
grade_concrete_label = Label(app, text="Grade of concrete", font=('normal', 10), padx=5, pady=5)
grade_concrete_label.grid(row=10, column=0)
grade_concrete_entry = Entry(app, textvariable=grade_concrete_text, width=30)
grade_concrete_entry.grid(row=10, column=1)
grade_concrete_unit = Label(app, text="N/mm2", font=('normal', 10), padx=5, pady=5);
grade_concrete_unit.grid(row=10, column=2);


grade_steel_text = IntVar()
grade_steel_text.set(0)
grade_steel_label = Label(app, text="Grade of longitudinal steel", font=('normal', 10), padx=5, pady=5)
grade_steel_label.grid(row=11, column=0)
grade_steel_entry = Entry(app, textvariable=grade_steel_text, width=30)
grade_steel_entry.grid(row=11, column=1)
grade_steel_unit = Label(app, text="N/mm2", font=('normal', 10), padx=5, pady=5);
grade_steel_unit.grid(row=11, column=2);


imposed_load_text = IntVar()
imposed_load_text.set(0)
imposed_load_label = Label(app, text="Imposed Load", font=('normal', 10), padx=5, pady=5)
imposed_load_label.grid(row=12, column=0)
imposed_load_entry = Entry(app, textvariable=imposed_load_text, width=30)
imposed_load_entry.grid(row=12, column=1)
imposed_load_unit = Label(app, text="kN", font=('normal', 10), padx=5, pady=5);
imposed_load_unit.grid(row=12, column=2);


no_legs_stirrups_text = IntVar()
no_legs_stirrups_text.set(0)
no_legs_stirrups_label = Label(app, text="Number of legs of Stirrups", font=('normal', 10), padx=5, pady=5)
no_legs_stirrups_label.grid(row=13, column=0)
no_legs_stirrups_entry = Entry(app, textvariable=no_legs_stirrups_text, width=30)
no_legs_stirrups_entry.grid(row=13, column=1)

shear_force_text = IntVar()
shear_force_text.set(0)
shear_force_label = Label(app, text="Shear Force", font=('normal', 10), padx=5, pady=5)
shear_force_label.grid(row=14, column=0)
shear_force_entry = Entry(app, textvariable=shear_force_text, width=30)
shear_force_entry.grid(row=14, column=1)
shear_force_unit = Label(app, text="kN", font=('normal', 10), padx=5, pady=5);
shear_force_unit.grid(row=14, column=2);

bending_moment_text = IntVar()
bending_moment_text.set(0)
bending_moment_label = Label(app, text="Bending Moment", font=('normal', 10), padx=5, pady=5)
bending_moment_label.grid(row=15, column=0)
bending_moment_entry = Entry(app, textvariable=bending_moment_text, width=30)
bending_moment_entry.grid(row=15, column=1)
bending_moment_unit = Label(app, text="kNm", font=('normal', 10), padx=5, pady=5);
bending_moment_unit.grid(row=15, column=2);
# bending_moment_label = Label(app, text="Non factoring BM", font=('normal', 10), padx=5, pady=5)
# bending_moment_label.grid(row=15, column=3)

elastic_modulus = 210000  # N/mm2
density_concrete = 25  # in KN/m^3


def save_var_latex(key, value):
    dict_var = {}

    file_path = os.path.join(os.getcwd(), "mydata.dat")

    try:
        with open(file_path, newline="") as file:
            reader = csv.reader(file)
            for row in reader:
                dict_var[row[0]] = row[1]
    except FileNotFoundError:
        pass

    dict_var[key] = value

    with open(file_path, "w") as f:
        for key in dict_var.keys():
            f.write(f"{key},{dict_var[key]}\n")

            
def pdfconverter():
    # TeX source filename
    tex_filename = 'BeamDesignReport.tex'
    filename, ext = os.path.splitext(tex_filename)
    # the corresponding PDF filename
    pdf_filename = filename + '.pdf'

    # compile TeX file
    subprocess.run(['pdflatex', '-interaction=nonstopmode', tex_filename])

    # check if PDF is successfully generated
    if not os.path.exists(pdf_filename):
        raise RuntimeError('PDF output not found')

    # open PDF with platform-specific command
    if platform.system().lower() == 'darwin':
        subprocess.run(['open', pdf_filename])
    elif platform.system().lower() == 'windows':
        os.startfile(pdf_filename)
    elif platform.system().lower() == 'linux':
        subprocess.run(['xdg-open', pdf_filename])
    else:
        raise RuntimeError('Unknown operating system "{}"'.format(platform.system()))


def getvals():
    # print(f"{clear_span_text.get()+width_Beam_text.get()}");
    effective_cover = clear_cover_text.get() + dia_tension_reinforcement_text.get() + dia_stirrups_text.get()
    effective_depth = overall_depth_text.get() - effective_cover
    effective_span = min(clear_span_text.get() + effective_depth / 1000,
                         clear_span_text.get() + (width_support_text.get() / 2000) + (width_support_text.get() / 2000))
    dead_load = (width_Beam_text.get() / 1000) * (overall_depth_text.get() / 1000) * density_concrete  # in KN/m
    factored_load = 1.5 * (dead_load + imposed_load_text.get())
    bending_moment = 1.5*bending_moment_text.get() + (factored_load * (effective_span ** 2)) / 8   # in KN-m
    shear_force = 1.5*shear_force_text.get() + factored_load * effective_span / 2

    def f(Xu_lim):
        return Xu_lim / (0.0035 * elastic_modulus) - (effective_depth - Xu_lim) / (
                0.002 * elastic_modulus + 0.87 * grade_steel_text.get())
        # equation for finding limiting value of Neutral axis depth

    X = 1.0  # initial guess
    Xu_lim = fsolve(f, X)  # solve the equation for Xu_lim

    lim_bending_moment = 0.36 * grade_concrete_text.get() * width_Beam_text.get() * Xu_lim * (effective_depth - (0.42 * Xu_lim)) / (10 ** 6)

    if bending_moment <= lim_bending_moment:
        print("Section should be designed as SIMPLY REINFORCED BEAM")

    else:
        print("Section should be designed as DOUBLY REINFORCED BEAM")

        def f(Ast_1):
            return Ast_1 - 0.414 * grade_concrete_text.get() / grade_steel_text.get() * width_Beam_text.get() * Xu_lim

        X = 1.0
        Ast_1 = fsolve(f, X)

        MR_2 = (bending_moment - lim_bending_moment)  # in KN-m

        def f(Ast_2):
            return MR_2 - 0.87 * grade_steel_text.get() * Ast_2 * (effective_depth - effective_cover) / 10 ** 6

        X = 1.0
        Ast_2 = fsolve(f, X)
        Ast = Ast_1 + Ast_2

        MAX_RF = 0.04 * width_Beam_text.get() * overall_depth_text.get()  # in mm^2
        MIN_RF = 0.85 * width_Beam_text.get() * effective_depth / grade_steel_text.get()

        if MAX_RF < Ast:
            print("Section is OVER REINFORCED , increase fy ")
        elif MIN_RF > Ast:
            Ast = MIN_RF
        else:
            print("THE REQUIRED REINFORCEMENT  ON TENSION SIDE IS ", Ast, " mm^2")

        n_TENSION_r = Ast / (3.14 * (dia_tension_reinforcement_text.get()/ 2) ** 2)  
        n_TENSION = int(Ast / (math.pi * (dia_tension_reinforcement_text.get()/ 2) ** 2)) + 1  # n= Number of BARS
        print("Number of BARS = ", n_TENSION, "of ", dia_tension_reinforcement_text.get(), "mm")
        Ast_provided = n_TENSION * math.pi * math.pow(dia_tension_reinforcement_text.get() / 2, 2)
        # print("Ast_provided", Ast_provided)
        pt_provided = Ast_provided * 100 / (width_Beam_text.get() * effective_depth)

        # print("pt_provided", pt_provided)

        def f(eps_sc):
            return eps_sc - (Xu_lim - effective_cover) * 0.0035 / Xu_lim

        X = 1.0
        eps_sc = fsolve(f, X)
        # print("eps_sc", eps_sc)

        if grade_steel_text.get() == 415:
            if eps_sc < 0.00144:
                fsc = eps_sc * elastic_modulus
            else:
                eps_sc_std = [0.00144, 0.00163, 0.00192, 0.00241, 0.002760, 0.00380]
                fsc_std = [288.7, 306.7, 324.8, 342.8, 351.8, 360.9]
                fsc = np.interp(eps_sc, eps_sc_std, fsc_std)

        else:
            if eps_sc < 0.00174:
                fsc = eps_sc * elastic_modulus
            else:
                eps_sc_std = [0.00174, 0.00195, 0.00226, 0.00277, 0.00312, 0.00417]
                fsc_std = [347.8, 369.6, 391.3, 413.0, 423.9, 434.8]
                fsc = np.interp(eps_sc, eps_sc_std, fsc_std)

        # print("fsc=", fsc)

        def f(Asc):
            return Asc * (fsc - 0.45 * grade_concrete_text.get()) - 0.87 * grade_steel_text.get() * Ast_2

        X = 1.0
        Asc = fsolve(f, X)
        # print("Asc=", Asc)

        print("THE REQUIRED REINFORCEMENT  ON COMPRESSION SIDE IS ", Asc, " mm^2")
        
        n_COMPRESSION_required = Asc / (math.pi * (dia_comp_reinforcement_text.get() / 2) ** 2)

        n_COMPRESSION = int(Asc / (math.pi * (dia_comp_reinforcement_text.get() / 2) ** 2)) + 1  # n= Number of BARS
        print("Number of BARS = ", n_COMPRESSION, "of ", dia_comp_reinforcement_text.get(), "mm")
        Asc_provided = n_COMPRESSION * math.pi * math.pow(dia_comp_reinforcement_text.get() / 2, 2)
        # print("Asc_provided=", Asc_provided)
        pc = Asc * 100 / (width_Beam_text.get() * effective_depth)
        pc_provided = Asc_provided * 100 / (width_Beam_text.get() * effective_depth)
        # print("pc_provided=", pc_provided)

    Design_Shear_Force = shear_force - (factored_load * effective_depth / 1000)
    # print("Vdesign = ", Design_Shear_Force)
    Nominal_Shear_Stress = (Design_Shear_Force * 1000) / (width_Beam_text.get() * effective_depth)
    # print("Tnominal= ", Nominal_Shear_Stress)

    fck_std = [10, 15, 20, 25, 30, 35, 40, 45, 50]
    Max_Shear_Capacity_std = [1.2, 2.0, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2]
    Max_Shear_Capacity = np.interp(grade_concrete_text.get(), fck_std, Max_Shear_Capacity_std)
    # print("Tcmax= ", Max_Shear_Capacity)

    if Nominal_Shear_Stress < Max_Shear_Capacity:
        fck_std = [15, 20, 25, 30, 35, 40]
        pt_std = [0.15, 0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00,4.00]
        Shear_Stress_std = np.array([(0.18, 0.18, 0.19, 0.20, 0.20, 0.20), (0.22, 0.22, 0.23, 0.23, 0.23, 0.23),
                                     (0.29, 0.30, 0.31, 0.31, 0.31, 0.32), (0.34, 0.35, 0.36, 0.37, 0.37, 0.38),
                                     (0.37, 0.39, 0.40, 0.41, 0.42, 0.42), (0.40, 0.42, 0.44, 0.45, 0.45, 0.46),
                                     (0.42, 0.45, 0.46, 0.48, 0.49, 0.49), (0.44, 0.47, 0.49, 0.50, 0.52, 0.52),
                                     (0.44, 0.49, 0.51, 0.53, 0.54, 0.55), (0.44, 0.51, 0.53, 0.55, 0.56, 0.57),
                                     (0.44, 0.51, 0.55, 0.57, 0.58, 0.60), (0.44, 0.51, 0.56, 0.58, 0.60, 0.62),
                                     (0.44, 0.51, 0.57, 0.60, 0.62, 0.63), (0.44, 0.51, 0.57, 0.60, 0.62, 0.63)])

        f = interp2d(fck_std, pt_std, Shear_Stress_std, kind='linear', fill_value=-1)
        Shear_Stress = f(grade_concrete_text.get(), pt_provided)
        # Shear_Stress = The permissible shear stress in concrete in beams without shear reinforcement
        SF = (Nominal_Shear_Stress - Shear_Stress)*effective_depth*width_Beam_text.get()/1000
        Asv = no_legs_stirrups_text.get() * (math.pi * (dia_stirrups_text.get() / 2) ** 2)
        Sv = (int((0.87 * grade_stirrups_text.get() * Asv * effective_depth / ((Nominal_Shear_Stress - Shear_Stress)*effective_depth *width_Beam_text.get())) / 10) + 1) * 10
        Asv_min = width_Beam_text.get() * Sv * (0.4 / (0.87 * grade_steel_text.get()))
        Spacing_Stirrup = (int((min(Sv, 300, 0.75 * effective_depth,
                                   0.87 * grade_stirrups_text.get() * Asv / (0.4 * width_Beam_text.get()))) / 10) - 1) * 10

        if Nominal_Shear_Stress < Shear_Stress:
            print("Nominal shear stirrups are provided  of ", Asv_min, " mm^2 Reinforcement")
        else:
            print("Section designed for shear force = ",
                  (Nominal_Shear_Stress - Shear_Stress) * width_Beam_text.get() * effective_depth / 1000,
                  "KN")
            print("Provide ", dia_stirrups_text.get(), "mm dia ", no_legs_stirrups_text.get(), "legged stirrups at ", Spacing_Stirrup, "mm spacing")

    else:
        print("Increase Depth of Beam")

    def getcheck():
        # Check for Tension Reinforcement
        check_Tension_RF_label = Label(app, text="Check for Tension Reinforcement", font=('normal', 10), padx=5,
                                       pady=5);
        # check_Tension_RF_label.pack()
        check_Tension_RF_label.grid(row=19, column=0);
        if pt_provided <= 4 and pt_provided > 0.85 * 100 / grade_steel_text.get():
            check_Tension_label = Label(app, text="Safe", font=('normal', 10), padx=5, pady=5);
            # check_Tension_label.pack()
            check_Tension_label.grid(row=19, column=1);
        else:
            check_Tension_label = Label(app, text="Unsafe", font=('normal', 10), padx=5, pady=5);
            # check_Tension_label.pack()
            check_Tension_label.grid(row=19, column=1);

        # Check for compression Reinforcement
        check_Compression_RF_label = Label(app, text="Check for Compression Reinforcement", font=('normal', 10), padx=5,
                                           pady=5);
        # check_Compression_RF_label.pack()
        check_Compression_RF_label.grid(row=20, column=0);
        if pc_provided <= 4 and pc_provided > 0.2:
            check_Compression_label = Label(app, text="Safe", font=('normal', 10), padx=5, pady=5);
            # check_Compression_label.pack()
            check_Compression_label.grid(row=20, column=1);
        else:
            check_Compression_label = Label(app, text="Unsafe", font=('normal', 10), padx=5, pady=5);
            # check_Compression_label.pack()
            check_Compression_label.grid(row=20, column=1);

        # Check for shear
        check_Shear_RF_label = Label(app, text="Check for Shear", font=('normal', 10), padx=5, pady=5);
        # check_Shear_RF_label.pack()
        check_Shear_RF_label.grid(row=21, column=0);
        if Nominal_Shear_Stress < Max_Shear_Capacity:
            check_Shear_label = Label(app, text="Safe", font=('normal', 10), padx=5, pady=5);
            # check_Shear_label.pack()
            check_Shear_label.grid(row=21, column=1);
        else:
            check_Shear_label = Label(app, text="Unsafe", font=('normal', 10), padx=5, pady=5);
            # check_Shear_label.pack()
            check_Shear_label.grid(row=21, column=1);

        def getoutput():
            clear_span_label.destroy()
            clear_span_entry.destroy()
            clear_span_unit.destroy()
            clear_cover_label.destroy()
            clear_cover_entry.destroy()
            clear_cover_unit.destroy()
            imposed_load_label.destroy()
            imposed_load_entry.destroy()
            imposed_load_unit.destroy()
            grade_steel_label.destroy()
            grade_steel_entry.destroy()
            grade_steel_unit.destroy()
            grade_concrete_label.destroy()
            grade_concrete_entry.destroy()
            grade_concrete_unit.destroy()
            grade_stirrups_label.destroy()
            grade_stirrups_entry.destroy()
            grade_stirrups_unit.destroy()
            no_legs_stirrups_label.destroy()
            no_legs_stirrups_entry.destroy()
            dia_stirrups_label.destroy()
            dia_stirrups_entry.destroy()
            dia_stirrups_unit.destroy()
            overall_depth_label.destroy()
            overall_depth_entry.destroy()
            overall_depth_unit.destroy()
            width_Beam_label.destroy()
            width_Beam_entry.destroy()
            width_Beam_unit.destroy()
            width_support_label.destroy()
            width_support_entry.destroy()
            width_support_unit.destroy()
            dia_tension_reinforcement_label.destroy()
            dia_tension_reinforcement_entry.destroy()
            dia_tension_reinforcement_unit.destroy()
            dia_comp_reinforcement_label.destroy()
            dia_comp_reinforcement_entry.destroy()
            dia_comp_reinforcement_unit.destroy()
            shear_force_label.destroy()
            shear_force_entry.destroy()
            shear_force_unit.destroy()
            bending_moment_label.destroy()
            bending_moment_entry.destroy()
            bending_moment_unit.destroy()
            button1.destroy()
            button2.destroy()
            check_Tension_RF_label.destroy()
            check_Compression_RF_label.destroy()
            check_Shear_label.destroy()
            check_Tension_label.destroy()
            check_Compression_label.destroy()
            check_Shear_RF_label.destroy()


            tab_control = ttk.Notebook(app, style="lefttab.TNotebook")

            tab1 = ttk.Frame(tab_control)
            tab2 = ttk.Frame(tab_control)
            tab3 = ttk.Frame(tab_control)
            tab4 = ttk.Frame(tab_control)
            tab5 = ttk.Frame(tab_control)
            tab6 = ttk.Frame(tab_control)

            # Add tabs to Notebook
            tab_control.add(tab1, text=f'{" Input Data ":^20s}')
            tab_control.add(tab2, text=f'{" Shear Force and Bending Moment ":^20s}')
            tab_control.add(tab3, text=f'{" Tension Reinforcement ":^20s}')
            tab_control.add(tab4, text=f'{" Compression Reinforcement ":^20s}')
            tab_control.add(tab5, text=f'{" Shear Design ":^20s}')
            tab_control.add(tab6, text=f'{" Export ":^20s}')

            # tab_control.pack(expand=1, fill='both')
            tab_control.grid(row=0,column=0, padx=5, pady=5)

            # label1 = Label(tab1, text='Input Data', font=('bold', 20), padx=20, pady=20)
            # label1.grid(column=0, row=0)
            #
            # label1 = Label(tab2, text='Shear Force and Bending Moment', font=('bold', 20), padx=20, pady=20)
            # label1.grid(column=0, row=0)
            #
            # label1 = Label(tab3, text='Tension Reinforcement', font=('bold', 20), padx=20, pady=20)
            # label1.grid(column=0, row=0)
            #
            # label1 = Label(tab4, text='Compression Reinforcement', font=('bold', 20), padx=20, pady=20)
            # label1.grid(column=0, row=0)
            #
            # label1 = Label(tab5, text='Compression Reinforcement', font=('bold', 20), padx=20, pady=20)
            # label1.grid(column=0, row=0)

            wrapper1 = LabelFrame(tab1, text="Input Data")
            wrapper2 = LabelFrame(tab2, text="Shear Force and Bending Moment")
            wrapper3 = LabelFrame(tab3, text="Tension Reinforcement")
            wrapper4 = LabelFrame(tab4, text="Compression Reinforcement")
            wrapper5 = LabelFrame(tab5, text="Shear Design")
            wrapper6 = LabelFrame(tab6, text="Generate Pdf")

            wrapper1.pack(fill="both", expand="YES", padx=20, pady=10)
            wrapper2.pack(fill="both", expand="YES", padx=20, pady=10)
            wrapper3.pack(fill="both", expand="YES", padx=20, pady=10)
            wrapper4.pack(fill="both", expand="YES", padx=20, pady=10)
            wrapper5.pack(fill="both", expand="YES", padx=20, pady=10)
            wrapper6.pack(fill="both", expand="YES", padx=20, pady=10)

            tree = ttk.Treeview(wrapper1, columns=(1, 2, 3, 4), show='headings', height='15')
            tree.pack()

            tree.heading(1, text="Data")
            tree.heading(2, text="Notation")
            tree.heading(3, text="Value")
            tree.heading(4, text="Unit")

            # Insert the data in Treeview widget
            tree.insert('', 'end', text="1", values=('Width of Beam', 'B', width_Beam_text.get(), 'mm'))
            if clear_span_text.get() != 0:
                tree.insert('', 'end', text="1", values=('Clear Span', 'L', clear_span_text.get(), 'm'))
            tree.insert('', 'end', text="1", values=('Width of Support', 'w', width_support_text.get(), 'mm'))
            tree.insert('', 'end', text="1", values=('Clear Cover', 'cc', clear_cover_text.get(), 'mm'))
            tree.insert('', 'end', text="1",
                        values=('Dia of Tension Reinforcement', 'dt', dia_tension_reinforcement_text.get(), 'mm'))
            tree.insert('', 'end', text="1",
                        values=('Dia of Compression Reinforcement', 'dc', dia_comp_reinforcement_text.get(), 'mm'))
            tree.insert('', 'end', text="1", values=('Dia of Stirrups', 'ds', dia_stirrups_text.get(), 'mm'))
            tree.insert('', 'end', text="1", values=('Grade of Concrete', 'fck', grade_concrete_text.get(), 'mm'))
            tree.insert('', 'end', text="1", values=('Grade of Steel', 'fy', grade_steel_text.get(), 'mm'))
            tree.insert('', 'end', text="1", values=('Grade of Stirrups', 'fys', grade_stirrups_text.get(), 'mm'))
            tree.insert('', 'end', text="1", values=('Effective Span', 'B', effective_span, 'm'))


            tree2 = ttk.Treeview(wrapper2, columns=(1, 2, 3, 4), show='headings', height='10')
            tree2.pack()

            tree2.heading(1, text="Data")
            tree2.heading(2, text="Notation")
            tree2.heading(3, text="Value")
            tree2.heading(4, text="Unit")

            tree2.insert('', 'end', text="1", values=('Self Wt', 'DL', float(f'{dead_load:.2f}'), 'kN/m'))
            tree2.insert('', 'end', text="1", values=('Imposed Load', 'LL', imposed_load_text.get(), 'kN/m'))
            tree2.insert('', 'end', text="1", values=('Factored Load', 'FL', float(f'{factored_load:.2f}'), 'kN/m'))
            tree2.insert('', 'end', text="1", values=('Factored Shear Force', 'Vu', float(f'{shear_force:.2f}'), 'kN'))
            tree2.insert('', 'end', text="1",
                         values=('Factored Bending Moment', 'Mu', float(f'{bending_moment:.2f}'), 'kN-m'))
            tree2.insert('', 'end', text="1",
                         values=('Limiting Bending Moment', 'Mu-lim', float(f'{lim_bending_moment[0]:.2f}'), 'kN-m'))

            #    tree2.grid(row=14, column=0, sticky='nsew')

            tree3 = ttk.Treeview(wrapper3, columns=(1, 2, 3, 4), show='headings', height='10')
            tree3.pack()

            tree3.heading(1, text="Data")
            tree3.heading(2, text="Notation")
            tree3.heading(3, text="Value")
            tree3.heading(4, text="Unit")

            tree3.insert('', 'end', text="1", values=(' ', 'Ast1', float(f'{Ast_1[0]:.2f}'), 'mm2'))
            tree3.insert('', 'end', text="1", values=(' ', 'Ast2', float(f'{Ast_2[0]:.2f}'), 'mm2'))
            tree3.insert('', 'end', text="1", values=('Total Reinforcement', 'Ast', float(f'{Ast[0]:.2f}'), 'mm2'))
            tree3.insert('', 'end', text="1", values=('No of Bars Provided', 'n', n_TENSION, ' '))
            #tree3.insert('', 'end', text="1",
            #             values=('Percentage Reinforcement Required', 'pt_required', float(f'{pt:.2f}'), '%'))
            tree3.insert('', 'end', text="1",
                         values=('Percentage Reinforcement Provided', 'pt_provided', float(f'{pt_provided:.2f}'), '%'))

            #   tree3.grid(row=14, column=0, sticky='nsew')


            tree4 = ttk.Treeview(wrapper4, columns=(1, 2, 3, 4), show='headings', height='10')
            tree4.pack()

            tree4.heading(1, text="Data")
            tree4.heading(2, text="Notation")
            tree4.heading(3, text="Value")
            tree4.heading(4, text="Unit")

            tree4.insert('', 'end', text="1",
                         values=('Stress in compression steel', 'fsc', float(f'{fsc[0]:.2f}'), 'N/mm2'))
            tree4.insert('', 'end', text="1", values=('Total Reinforcement', 'Asc', float(f'{Asc[0]:.2f}'), 'mm2'))
            tree4.insert('', 'end', text="1", values=('No of Bars Provided', 'n', n_COMPRESSION, ' '))
            tree4.insert('', 'end', text="1",
                         values=('Percentage Reinforcement Required', 'pc_required', float(f'{pc[0]:.2f}'), '%'))
            tree4.insert('', 'end', text="1",
                         values=('Percentage Reinforcement Provided', 'pc_provided', float(f'{pc_provided:.2f}'), '%'))
            tree4.insert('', 'end', text="1",
                         values=('Reinforcement Provided', 'Asc_provided', float(f'{Asc_provided:.2f}'), 'mm2'))

            tree5 = ttk.Treeview(wrapper5, columns=(1, 2, 3, 4), show='headings', height='10')
            tree5.pack()

            tree5.heading(1, text="Data")
            tree5.heading(2, text="Notation")
            tree5.heading(3, text="Value")
            tree5.heading(4, text="Unit")

            tree5.insert('', 'end', text="1", values=('Factored Shear Force', 'Vu', float(f'{shear_force:.2f}'), 'kN'))
            tree5.insert('', 'end', text="1", values=('Design Shear Force', 'V_design', float(f'{Design_Shear_Force:.2f}'), 'kN'))
            tree5.insert('', 'end', text="1", values=('Nominal Shear Stress', 'τv', float(f'{Nominal_Shear_Stress:.2f}'), 'N/mm2'))
            tree5.insert('', 'end', text="1", values=('Max Shear Stress', 'τcmax',float(f'{Max_Shear_Capacity:.2f}'), 'N/mm2'))
            tree5.insert('', 'end', text="1", values=('Area of Stirrups', 'Asv', float(f'{Asv:.2f}'), 'mm2'))
            tree5.insert('', 'end', text="1", values=('Spacing of Stirrups', 'Sv', float(f'{Spacing_Stirrup:.2f}'), 'mm'))

            def save_var_latex(key, value):
                dict_var = {}

                file_path = os.path.join(os.getcwd(), "mydata.dat")

                try:
                    with open(file_path, newline="") as file:
                        reader = csv.reader(file)
                        for row in reader:
                            dict_var[row[0]] = row[1]
                except FileNotFoundError:
                    pass

                dict_var[key] = value

                with open(file_path, "w") as f:
                    for key in dict_var.keys():
                        f.write(f"{key},{dict_var[key]}\n")
            
            button3.destroy()
            save_var_latex("clear_span", clear_span_text.get())
            save_var_latex("width_Beam", width_Beam_text.get())
            save_var_latex("width_support", width_support_text.get())
            save_var_latex("cover", clear_cover_text.get())
            save_var_latex("dia_tension_reinforcement", dia_tension_reinforcement_text.get())
            save_var_latex("dia_comp_reinforcement", dia_comp_reinforcement_text.get())
            save_var_latex("dia_stirrups", dia_stirrups_text.get())
            save_var_latex("grade_concrete", grade_concrete_text.get())
            save_var_latex("grade_steel", grade_steel_text.get())
            save_var_latex("grade_stirrups", grade_stirrups_text.get())
            save_var_latex("effective_span", effective_span)
            save_var_latex("effective_cover", effective_cover)
            save_var_latex("effective_depth", effective_depth)
            save_var_latex("dead_load", float(f'{dead_load:.2f}'))
            save_var_latex("imposed_load", imposed_load_text.get())
            save_var_latex("factored_load", float(f'{factored_load:.2f}'))
            save_var_latex("shear_force", float(f'{shear_force:.2f}'))
            save_var_latex("bending_moment", float(f'{bending_moment:.2f}'))
            save_var_latex("lim_bending_moment", float(f'{lim_bending_moment[0]:.2f}'))
            save_var_latex("Ast_1", float(f'{Ast_1[0]:.2f}'))
            save_var_latex("Ast_2", float(f'{Ast_2[0]:.2f}'))
            save_var_latex("Ast", float(f'{Ast[0]:.2f}'))
            save_var_latex("n_TENSION", n_TENSION)
            save_var_latex("pt_provided", float(f'{pt_provided:.2f}'))
            save_var_latex("fsc", float(f'{fsc[0]:.2f}'))
            save_var_latex("Asc", float(f'{Asc[0]:.2f}'))
            save_var_latex("n_COMPRESSION", n_COMPRESSION)
            save_var_latex("pc", float(f'{pc[0]:.2f}'))
            save_var_latex("pc_provided", float(f'{pc_provided:.2f}'))
            save_var_latex("Asc_provided", float(f'{Asc_provided:.2f}'))
            save_var_latex("Vu", float(f'{shear_force:.2f}'))
            save_var_latex("Design_Shear_Force", float(f'{Design_Shear_Force:.2f}'))
            save_var_latex("Nominal_Shear_Stress", float(f'{Nominal_Shear_Stress:.2f}'))
            save_var_latex("Max_Shear_Capacity", float(f'{Max_Shear_Capacity:.2f}'))
            save_var_latex("Asv", float(f'{Asv:.2f}'))
            save_var_latex("Sv", float(f'{Spacing_Stirrup:.2f}'))
            save_var_latex("n_TENSION_required", float(f'{n_TENSION_r[0]:.2f}'))
            save_var_latex("overall_depth", overall_depth_text.get())
            save_var_latex("MIN_RF", float(f'{MIN_RF:.2f}'))
            save_var_latex("MAX_RF", float(f'{MAX_RF:.2f}'))
            save_var_latex("Xu_lim", float(f'{Xu_lim[0]:.2f}'))
            save_var_latex("esc", float(f'{eps_sc[0]:.5f}'))
            save_var_latex("n_COMPRESSION_required", float(f'{n_COMPRESSION_required[0]:.2f}'))
            save_var_latex("Xu_lim", float(f'{Xu_lim[0]:.2f}'))
            save_var_latex("Shear_Stress", float(f'{Shear_Stress[0]:.2f}'))
            save_var_latex("SF", float(f'{SF[0]:.2f}'))
            save_var_latex("no_legs_stirrups", float(f'{no_legs_stirrups_text.get():.2f}'))
            save_var_latex("Sv_1", float(f'{Sv:.2f}'))
            
            
            def pdfconverter():
                with open('BeamDesignReport.tex','w') as file:
                    file.write('\\documentclass[12pt,a4paper]{article}\n')
                    file.write('\\usepackage{marginnote} \n')
                    file.write('\\usepackage{wallpaper}\n')
                    file.write('\\usepackage{lastpage}\n')
                    file.write('\\usepackage[left=1.3cm,right=4.6cm,top=1.8cm,bottom=4.0cm,marginparwidth=3.4cm]{geometry}\n ')
                    file.write('\\usepackage{amsmath}\n')
                    file.write('\\usepackage{amssymb}\n')
                    file.write('\\usepackage{xcolor}\n')
                    file.write('\\usepackage{datatool}\n')   
                    file.write('\\usepackage{filecontents}\n')   
                    file.write('\\DTLsetseparator{,}\n')   
                    file.write('\\DTLloaddb[nosheader, keys={thekey,thevalue}]{mydata}{mydata.dat}\n') 
                    file.write('\\newcommand{\\var}[1]{\DTLfetch{mydata}{thekey}{#1}{thevalue}}\n') 
                    file.write('\\usepackage{fancyhdr}\n') 
                    file.write('\\setlength{\headheight}{80pt} \n') 
                    file.write('\\pagestyle{fancy}\\fancyhf{}\n ') 
                    file.write('\\renewcommand{\headrulewidth}{0pt} \n') 
                    file.write('\\setlength{\parindent}{0cm}\n') 
                    file.write('\\newcommand{\\tab}{\hspace*{2em}}\n ') 
                    file.write('\\newcommand\BackgroundStructure{ \n') 
                    file.write('\\setlength{\\unitlength}{1mm}\n') 
                    file.write(' \\setlength\\fboxsep{0mm} \n') 
                    file.write('\\setlength\\fboxrule{0.5mm}\n') 
                    file.write('\put(10, 10){\\fcolorbox{black}{blue!10}{\\framebox(155,247){}}}\n') 
                    file.write('\put(165, 10){\\fcolorbox{black}{blue!10}{\\framebox(37,247){}}}\n') 
                    file.write('\put(10, 262){\\fcolorbox{black}{white!10}{\\framebox(192, 25){}}}\n') 
                    file.write('\put(137, 263){\includegraphics[height=23mm,keepaspectratio]{logo}}\n') 
                    file.write('}\n') 
                    file.write('\\fancyhead[L]{\\begin{tabular}{l r | l r}\n') 
                    file.write('\\textbf{Project} & BEAM DESIGN & \\textbf{Page} & \\thepage/\pageref{LastPage} \\\ \n') 
                    file.write('\\textbf{Designer} & Saurabh k. Khandelwal & \\textbf{Reviewer} & PROF. M. N. Shariff \\\ \n') 
                    file.write('\\textbf{Updated} & 27/11/2012 & \\textbf{Reviewed} & 2/12/2012 \\\ \n') 
                    file.write('\end{tabular}}\n') 
                    file.write('\\begin{document}\n') 
                    file.write('\AddToShipoutPicture{\BackgroundStructure} \n') 
                    file.write('\section{Design of Doubly Reinforced Beam} \n') 
                    file.write('\subsection{Preliminary Design} \n\n\n') 
                    file.write('Sectional characteristic:\n\n') 
                    file.write('\\tab Width\; of\; the\; beam ($b$)\; =\; \\var{width_Beam}mm \\\[8pt]\n')
                    file.write('\\tab Overall\; depth\; of\; the\; beam ($D$)\; =\; \\var{overall_depth}mm \\\[8pt]\n')
                    file.write('\\tab Width\; of\; the\; support ($w$)\; =\; \\var{width_support}mm \\\[8pt]\n')
                    file.write('\\tab Clear\; cover\; =\; \\var{cover} mm \\\[8pt]\n')
                    file.write('\\tab Diameter\; of\; the\; tension\; bar ($d_t$)\; =\; \\var{dia_tension_reinforcement}mm \\\[8pt]\n')
                    file.write('\\tab Diameter\; of\; the\; compression\; bar ($d_c$)\; =\; \\var{dia_comp_reinforcement}mm \\\[8pt]\n')
                    file.write('\\tab Diameter\; of\; the\; stirrup\; bar ($d_s$)\; =\; \\var{dia_stirrups}mm \\\[8pt]\n')   
                    file.write('\\tab Effective\; depth(d)\; =\; D-clear\; cover-\\frac{d_t}{2}-d_s \\\[8pt]\n\n') 
                    file.write('\\begin{equation*}\n') 
                    file.write('d = \\var{overall_depth}-\\var{cover}-\\frac{\\var{dia_tension_reinforcement}}{2}-\\var{dia_stirrups}\n') 
                    file.write('\end{equation*}\n\n') 
                    file.write('\\begin{equation*}\n') 
                    file.write('d = \\var{effective_depth}mm\n') 
                    file.write('\end{equation*}\n\n')
#                     file.write('\\tab d\; =\; \\var{overall_depth}-\\var{clear_cover}-\\frac{\\var{dia_tension_reinforcement}}{2}-\\var{dia_stirrups} \\\[8pt]\n')
#                     file.write('\\tab d\; =\; \\var{effective_depth}mm \\\[8pt]\n')
                    if clear_span_text.get() != 0 :
                        file.write('\\tab Clear\; span(l)\; =\; \\var{clear_span}m\\\[8pt]\n') 
                        file.write('\\tab Effective\; span({l_e}_f_f)\; =\; min \\begin{Bmatrix}l+d = \\var{clear_span} + \\var{effective_depth}  \\\l+w = \\var{clear_span} + \\var{width_support} \end{Bmatrix} = \\var{effective_span} m \\\[8pt]\n\n')
                    file.write('Material characteristic:\n\n\n') 
                    file.write('\\tab Grade\; of\; concrete({f_c}_k)\; =\; \\var{grade_concrete}\\frac{N}{mm^2} \\\[8pt]\n')
                    file.write('\\tab Grade\; of\; steel(f_y)\; =\; \\var{grade_steel}\\frac{N}{mm^2} \\\[8pt]\n')
                    file.write('\\tab Grade\; of\; stirrup(f_s)\; =\; \\var{grade_stirrups}\\frac{N}{mm^2} \\\[8pt]\n\n')
                    if clear_span_text.get() != 0 :
                        file.write('Grading characteristic:\n\n') 
                        file.write('\\begin{equation*}\n Dead\; load(DL) = \\var{dead_load}\\frac{kN}{m}\n \end{equation*}\n\n')
                        file.write('\\begin{equation*}\n Imposed\; load(LL) = \\var{imposed_load}\\frac{kN}{m}\n \end{equation*}\n\n')
                        file.write('\\begin{equation*}\n Total\; Factored\; load(w_u) = 1.5\\times \left ( DL + LL \\right ) = \\var{factored_load}\\frac{kN}{m}\n \end{equation*}\n\n')
#                         file.write('\\tab Dead\; load(DL)\; =\; \\var{dead_load}\\frac{kN}{m} \\\[8pt]\n')
#                         file.write('\\tab Imposed\; load(LL)\; =\; \\var{imposed_load}\\frac{kN}{m} \\\[8pt]\n')
#                         file.write('\\tab Total\; Factored\; load(w_u)\; =\; 1.5\\times \left ( DL + LL \\right )\; =\; \\var{factored_load}\\frac{kN}{m} \\\[8pt]\n\n')
                    file.write('\subsection{Preliminary Design} \n\n') 
                    file.write('Shear Force($V_u$):\n\n') 
                    if clear_span_text.get() != 0 :
                        file.write('\\begin{equation*}\n') 
                        file.write('V_u = \\frac{w_u{l_e}_f_f}{2}\n') 
                        file.write('\end{equation*}\n\n')
                        file.write('\\begin{equation*}\n') 
                        file.write('V_u = \\frac{\\var{factored_load}\\times \\var{effective_span}}{2}\n') 
                        file.write('\end{equation*}\n\n')
                    file.write('\\begin{equation*}\n') 
                    file.write('V_u = \\var{shear_force} kN\n') 
                    file.write('\end{equation*}\n\n')    
                    file.write('Ultimate bending moment:\n') 
                    if clear_span_text.get() != 0 :
                        file.write('\\begin{equation*}\n') 
                        file.write('M_u = \\frac{w_u{l_e}_f_f^2}{8}\n') 
                        file.write('\end{equation*}\n\n') 
                        file.write('\\begin{equation*}\n') 
                        file.write('M_u = \\frac{\\var{factored_load}\\times \\var{effective_span}^2}{8}\n') 
                        file.write('\end{equation*}\n\n') 
                    file.write('\\begin{equation*}\n') 
                    file.write('M_u = \\var{bending_moment} kNm\n') 
                    file.write('\end{equation*}\n\n') 
                    file.write('Limiting bending moment:\n\n') 
                    file.write('\\tab From\; similiar\; triangle\\\[8pt]\n')
                    file.write('\\begin{equation*}\n')
                    file.write('\\frac{{x_u}_l_i_m}{0.0035} = \\frac{d-{x_{u,lim}}}{0.002+\\frac{0.87f_y}{E_s}}\n')
                    file.write('\end{equation*}\n')
                    file.write('\\begin{equation*} \n{x_u}_l_i_m = \\var{Xu_lim}mm \n\end{equation*}\n')
                    file.write('\\begin{equation*} \n {M_u}_l_i_m = 0.36f_c_kb{x_u}_l_i_m\left ( d -0.42x_u\\right)\n \end{equation*}')
                    file.write('\\begin{equation*}\n')
                    file.write('{M_u}_l_i_m = 0.36\\times \\var{grade_concrete}\\times \\var{width_Beam}\\times \\var{Xu_lim}\\times \left (\\var{effective_depth} -0.42\\times \\var{Xu_lim}\\right )\n')
                    file.write('\end{equation*}\n\n')
                    file.write('\\begin{equation*}\n {M_u}_l_i_m = \\var{lim_bending_moment}kNm < M_u\n \end{equation*}\n')
                    file.write('So, Section should be designed as Doubly reinforced beam.\n')
                    file.write('\par\\vspace{\\baselineskip}\n\n')
                    file.write('\subsubsection{Tension Reinforcement}\n\n')
                    file.write('Total tensile reinforcement:\n')
                    file.write('\\begin{equation*}\n A_s_t={A_s_t}_l_i_m+\Delta A_s_t \n \end{equation*}\n')
                    file.write('\\tab Where\; {A_s_t}_l_i_m\; is\; due\; to\; tensile\; steel\; corresponding\; to\; {M_u}_l_i_m \\\[8pt]\n')
                    file.write('\\tab and\; \Delta A_s_t\; is\; corresponding\; to\; \Delta M_u\\\[8pt]\n\n')
                    file.write('Limiting Tensile Reinforcement:\n')
                    file.write('\\begin{equation*}\n')
                    file.write('{A_s_t}_l_i_m = 0.5\\frac{f_c_k}{f_y}bd\left [ 1-\sqrt{1-\\frac{4.6{M_u}_l_i_m}{f_c_kbd^2}} \\right ]\n')
                    file.write('\end{equation*}\n')
                    file.write('\\begin{equation*}\n')
                    file.write('{A_s_t}_l_i_m = 0.5\\times \\frac{\\var{grade_concrete}}{\\var{grade_steel}}\\times \\var{width_Beam}\\times \\var{effective_depth}\\times \left [ 1-\sqrt{1-\\frac{4.6\\times \\var{lim_bending_moment}}{\\var{grade_concrete}\\times \\var{width_Beam}\\times \\var{effective_depth}^2}} \\right ]\n')
                    file.write('\end{equation*}\n') 
                    file.write('\\begin{equation*}\n {A_{st,lim}}= \\var{Ast_1}mm^2\n \end{equation*}\n')
                    file.write('\Delta M_u:\n\n')
                    file.write("\\begin{equation*}\n \Delta M_u = 0.87f_y\Delta A_s_t\left ( d-d{}' \\right )\n \end{equation*}")
                    file.write('\\begin{equation*}\n')
                    file.write('\\var{bending_moment}-\\var{lim_bending_moment} = 0.87\\times \\var{grade_steel}\\times \Delta A_s_t\\times \left (\\var{effective_depth}-\\var{effective_cover} \\right )\n')
                    file.write('\end{equation*}\n')
                    file.write('\\begin{equation*}\n \Delta A_s_t = \\var{Ast_2}mm^2\n \end{equation*}\n')
                    file.write("\\tab Where d' is effective cover\\\[8pt]")
                    file.write('\\begin{equation*}\n A_s_t=\\var{Ast_1}+\\var{Ast_2}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n A_s_t=\\var{Ast}mm^2\n \end{equation*}\n')
                    file.write('Number of Bars:\n\n')
                    file.write('\\begin{equation*}\n n = \\frac{A_s_t}{\\frac{\pi}{4}d_t^2}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n n = \\frac{\\var{Ast}}{\\frac{\pi}{4}\\times \\var{dia_tension_reinforcement}^2}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n n = \\var{n_TENSION_required}\n \end{equation*}\n')
                    file.write('\\tab So\; provide\; \\var{n_TENSION}\; bars\; of\; \\var{dia_tension_reinforcement}\; mm\; dia\; as\; tensile\; reinforcement.\; \\\[8pt]')
                    file.write('Percentage reinforcement:\n \n')
                    file.write('\\begin{equation*}\n pt = \\frac{n\\times \\frac{\pi}{4}d_t^2}{bd} \end{equation*}\n')
                    file.write('\\begin{equation*}\n pt = \\frac{\\var{n_TENSION}\\times \\frac{\pi}{4}\\times \\var{dia_tension_reinforcement}^2}{\\var{width_Beam}\\times \\var{effective_depth}}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n pt = \\var{pt_provided}\%\n \end{equation*}\n')
                    file.write('Check:\n\n')
                    file.write('\\begin{equation*}\n {A_s}_t_m_a_x = 0.04bD\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n {A_s}_t_m_a_x = 0.04\\times \\var{width_Beam}\\times \\var{overall_depth}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n {A_s}_t_m_a_x = \\var{MAX_RF} > A_s_t \n \end{equation*}\n \marginnote{OK}\n')
                    file.write('\\begin{equation*}\n {A_s}_t_m_i_n = \\frac{0.85bd}{f_y}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n {A_s}_t_m_i_n = \\frac{0.85\\times \\var{width_Beam}\\times \\var{effective_depth}}{\\var{grade_steel}}\n \end{equation*}\n') 
                    file.write('\\begin{equation*}\n {A_s}_t_m_i_n = \\var{MIN_RF} < A_s_t\n \end{equation*}\n \marginnote{OK}\n')
                    file.write('\\par\\vspace{\\baselineskip}\n')
                    file.write('\subsubsection{Compression Reinforcement}\n\n')
                    file.write("\\begin{equation*}\n \epsilon _s_c = \left(\\frac{{x_u}_l_i_m-d{}'}{{x_u}_l_i_m} \\right )\\times 0.0035\n \end{equation*}\n")
                    file.write('\\begin{equation*}\n \epsilon _s_c = \left (\\frac{\\var{Xu_lim}-\\var{effective_cover}}{Xu_lim} \\right )\\times 0.0035\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n \epsilon _s_c = \\var{esc}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n f_s_c = \\var{fsc}\\frac{N}{mm^2}\n \marginnote{(SP-16 Table A)}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n \left ( f_s_c - 0.45f_c_k \\right )A_s_c = 0.87f_y\Delta A_s_t\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n \left (\\var{fsc} - 0.45\\var{grade_concrete} \\right )\\times A_s_c = 0.87\\times \\var{grade_steel}\\times \\var{Ast_2}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n A_s_c = \\var{Asc}mm^2\n \end{equation*}\n')
                    file.write('Number of Bars:\n\n')
                    file.write('\\begin{equation*}\n n = \\frac{A_s_c}{\\frac{\pi}{4}d_c^2}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n n = \\frac{\\var{Asc}}{\\frac{\pi}{4}\\times \\var{dia_comp_reinforcement}^2}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n n = \\var{n_COMPRESSION_required}\n \end{equation*}\n')
                    file.write('\\tab So\; provide\; \\var{n_COMPRESSION}\; bars\; of\; \\var{dia_comp_reinforcement}\; mm\; dia\; as\; compressive\; reinforcement.\; \\\[8pt]\n')
                    file.write('Percentage reinforcement:\n\n')
                    file.write('\\begin{equation*}\n pc = \\frac{n\\times \\frac{\pi}{4}d_c^2}{bd}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n pc = \\frac{\\var{n_COMPRESSION}\\times \\frac{\pi}{4}\\times \\var{dia_comp_reinforcement}^2}{\\var{width_Beam}\\times \\var{effective_depth}}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n pc = \\var{pc_provided}\%\n \end{equation*}\n')
                    file.write('Check:\n\n')
                    file.write('\\begin{equation*}\n {A_s}_c_m_a_x = 0.04bD\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n {A_s}_c_m_a_x = 0.04\\times \\var{width_Beam}\\times \\var{overall_depth}\n \end{equation*}\n')
                    file.write('\\begin{equation*}\n {A_s}_c_m_a_x = \\var{MAX_RF} > A_s_c\n \end{equation*}\n \marginnote{OK}\n')
                    file.write('\subsubsection{Shear Reinforcement}\n\n') 
                    file.write('\\tab Design\; shear\; force(V_d)\; =\; \\var{Design_Shear_Force}kN \\\[8pt]\n')
                    file.write('\\tab Nominal\; shear\; stress(\\tau _v)\; =\; \\frac{V_d}{bd} \\\[8pt]\n\n')
                    file.write('\\begin{equation*}\n \\tau _v = \\frac{\\var{Design_Shear_Force}\\times 1000}{\\var{width_Beam}\\times \\var{effective_depth}}\n \end{equation*}\n\n')
                    file.write('\\begin{equation*}\n \\tau _v = \\var{Nominal_Shear_Stress}\\frac{N}{mm^2}\n \end{equation*}\n\n')
#                     file.write('\\tab (\\tau _v)\; =\; \\frac{\\var{Design_Shear_Force}\\times 1000}{\\var{width_Beam}\\times \\var{effective_depth}} \\\[8pt]\n')
#                     file.write('\\tab (\\tau _v)\; =\; \\var{Nominal_Shear_Stress}\\frac{N}{mm^2} \\\[8pt]\n')
                    file.write('\\tab Max\; shear\; capacity({\\tau _m}_a_x)\; =\; \\var{Max_Shear_Capacity}\\frac{N}{mm^2}\n \marginnote{OK} \\\[8pt]\n\n')
#                     file.write('\n')
#                     file.write('\marginnote{IS 456 : Table 20}\n')
                    file.write('\\tab Permissible\; shear\; capacity(\\tau _c)\; =\; \\var{Shear_Stress}\\frac{N}{mm^2} \marginnote{IS 456 : Table 19} \\\[8pt]\n')
#                     file.write('\marginnote{IS 456 : Table 19}\n')
                    if Nominal_Shear_Stress < Shear_Stress:
                        file.write('\\tab As\; \\tau _v <\\tau _c \\\[8pt]\n')
                        file.write('\\tab So\; Nominal\; shear\; stirrups\; are\; provided\; \\\[8pt]\n')
#                         file.write('')
#                         file.write('')
                    else :     
                        file.write('\\tab As\; \\tau _v >\\tau _c \\\[8pt]\n')
                        file.write('\\tab So\; shear\; reinforcement\; is\; designed\; for\; \\\[8pt]\n\n')
                        file.write('\\begin{equation*}\n Shear force(V_u_s) = \left ( \\tau _v -\\tau _c \\right )bd\n \end{equation*}\n')
                        file.write('\\begin{equation*}\n V_u_s = \left ( \\var{Nominal_Shear_Stress} -\\var{Shear_Stress} \\right )\\times \\var{width_Beam}\\times \\var{effective_depth}\\times \\frac{1}{1000}\n \end{equation*}\n')
                        file.write('\\begin{equation*}\n V_u_s = \\var{SF}kN\n \end{equation*}\n\n')
                        file.write('\\begin{equation*}\n A_s_v = n\\times \\frac{\pi }{4}d_s^2\n \end{equation*}\n')
                        file.write('\\tab here\; n\; is\; number\; of\; legs\; of\; stirrups \\\[8pt]\n')
                        file.write('\\begin{equation*}\n A_s_v = \\var{no_legs_stirrups}\\times \\frac{\pi }{4}\\var{dia_stirrups}^2\n \end{equation*}\n')
                        file.write('\\begin{equation*}\n A_s_v = \\var{Asv}mm^2\n \end{equation*}\n\n')
                        file.write('\\tab Spacing\; of\; stirrups\; \\\[8pt]\n')
                        file.write('\\begin{equation*}\n S_v = \\frac{0.87f_yA_s_vd}{V_u_s}\n \end{equation*}\n')
                        file.write('\\begin{equation*}\n S_v = \\frac{0.87\\times \\var{grade_stirrups}\\times \\var{Asv}\\times \\var{effective_depth}}{\\var{SF}}\n \end{equation*}\n')
                        file.write('\\begin{equation*}\n S_v = \\var{Sv_1}\n \end{equation*}\n')
                        file.write('\\begin{equation*}\n Spacing \leq min\\begin{Bmatrix}S_v = \\var{Sv_1}\\\ 300\; mm\\\ 0.75d = 0.75\\times \\var{effective_depth} \\\ \\frac{A_s_v}{bS_v}\geq \\frac{0.4}{0.87f_y} = \\frac{\\var{Asv}}{\\var{width_Beam}\\times Sv}\geq \\frac{0.4}{0.87\\times \\var{grade_steel}} \end{Bmatrix}\n \end{equation*}\n')
                        file.write('\\begin{equation*}\n Spacing = \\var{Sv}mm\n \end{equation*}\n\n')
                        file.write('\\tab So\; provide\; \\var{no_legs_stirrups}-legged\; \\var{dia_stirrups}mm\; stirrups\; @\; \\var{Sv}mm\; c/c  \\\[8pt]\n')
                    file.write('\end{document}\n')

                # TeX source filename
                tex_filename = 'BeamDesignReport.tex'
                filename, ext = os.path.splitext(tex_filename)
                # the corresponding PDF filename
                pdf_filename = filename + '.pdf'

                # compile TeX file
                subprocess.run(['pdflatex', '-interaction=nonstopmode', tex_filename])

                # check if PDF is successfully generated
                if not os.path.exists(pdf_filename):
                    raise RuntimeError('PDF output not found')

                # open PDF with platform-specific command
                if platform.system().lower() == 'darwin':
                    subprocess.run(['open', pdf_filename])
                elif platform.system().lower() == 'windows':
                    os.startfile(pdf_filename)
                elif platform.system().lower() == 'linux':
                    subprocess.run(['xdg-open', pdf_filename])
                else:
                    raise RuntimeError('Unknown operating system "{}"'.format(platform.system()))
                    
            
            button_gen_pdf = Button(wrapper6,text="Generate Pdf",width=12,bg='#03A9F4',fg='#fff',command= pdfconverter)
            button_gen_pdf.pack()


#                 data = [ '\n','\n','\n',
#                          '1) Data : \n'
#                          'Width of Beam ', '(B) = ', width_Beam_text.get(), ' mm','\n',
#                          'Clear Span ', '(L) = ', clear_span_text.get(), ' m','\n',
#                          'Width of Support ', ' (w) = ', width_support_text.get(), 'mm','\n',
#                          'Clear Cover = ', clear_cover_text.get(), ' mm','\n',
#                          'Dia of Tension Reinforcement = ', dia_tension_reinforcement_text.get(), ' mm','\n',
#                          'Dia of Compression Reinforcement = ', dia_comp_reinforcement_text.get(), ' mm','\n',
#                          'Dia of Stirrups = ', dia_stirrups_text.get(), ' mm','\n',
#                          'Grade of Concrete ', '(fck) = ', grade_concrete_text.get(), 'N/mm2','\n',
#                          'Grade of Steel ', '(fy) = ', grade_steel_text.get(), ' N/mm2','\n',
#                          'Grade of Stirrups = ', grade_stirrups_text.get(), ' N/mm2','\n',
#                          'Effective Span ', '(Leff) = ', effective_span, ' m','\n','\n',
#                          '2) Shear Force and Bending Moment:','\n',
#                          'Self Wt ', '(DL) = ', float(f'{dead_load:.2f}'), ' kN/m','\n',
#                          'Imposed Load ', '(LL) = ', imposed_load_text.get(), ' kN/m','\n',
#                          'Factored Load', ' (FL) = ', float(f'{factored_load:.2f}'), ' kN/m','\n',
#                          'Factored Shear Force', ' (Vu) = ', float(f'{shear_force:.2f}'), ' kN','\n',
#                          'Factored Bending Moment', ' (Mu) = ', float(f'{bending_moment:.2f}'), ' kN-m','\n',
#                          'Limiting Bending Moment', ' (Mu-lim) = ', float(f'{lim_bending_moment[0]:.2f}'), ' kN-m','\n','\n',
#                          '3) Tension Reinforcement:','\n',
#                          'Ast1 = ', float(f'{Ast_1[0]:.2f}'), ' mm2','\n',
#                          'Ast2 = ', float(f'{Ast_2[0]:.2f}'), ' mm2','\n',
#                          'Total Reinforcement', ' (Ast) = ', float(f'{Ast[0]:.2f}'), ' mm2','\n',
#                          'No of Bars Provided', ' (n) = ', n_TENSION, '\n',
#                          #'Percentage Reinforcement Required = ', float(f'{pt:.2f}'), ' %','\n',
#                          'Percentage Reinforcement Provided = ',float(f'{pt_provided:.2f}'), ' %','\n','\n',
#                          '\n','4) Compression Reinforcement:','\n',
#                          'Stress in compression steel', ' (fsc) = ', float(f'{fsc[0]:.2f}'), ' N/mm2','\n',
#                          'Total Reinforcement', ' (Asc) = ', float(f'{Asc[0]:.2f}'), ' mm2','\n',
#                          'No of Bars Provided', ' (n) = ', n_COMPRESSION,'\n',
#                          'Percentage Reinforcement Required = ', float(f'{pc[0]:.2f}'), ' %','\n',
#                          'Percentage Reinforcement Provided = ', float(f'{pc_provided:.2f}'), ' %','\n',
#                          'Reinforcement Provided = ', float(f'{Asc_provided:.2f}'), ' mm2','\n','\n','\n',
#                          '5) Shear Design:','\n',
#                          'Factored Shear Force', ' (Vu) = ', float(f'{shear_force:.2f}'), ' kN','\n',
#                          'Design Shear Force = ', float(f'{Design_Shear_Force:.2f}'), ' kN','\n',
#                          'Nominal Shear Stress = ',  float(f'{Nominal_Shear_Stress:.2f}'), ' N/mm2','\n',
#                          'Max Shear Stress = ',float(f'{Max_Shear_Capacity:.2f}'), ' N/mm2','\n',
#                          'Area of Stirrups', ' (Asv) = ', float(f'{Asv:.2f}'), ' mm2','\n',
#                          'Spacing of Stirrups', ' (Sv) = ', float(f'{Spacing_Stirrup:.2f}'), ' mm','\n',
#                          # '(τv) = ',, ' (τcmax) = '
#                         ]

#                 pdf = fpdf.FPDF(format='letter')
#                 pdf.add_page()
#                 pdf.set_font("Arial", size=18)
#                 pdf.write(5, '                                     BEAM DESIGN REPORT        ')
#                 pdf.set_font("Arial", size=12)
#                 for i in data:
#                     pdf.write(5, str(i))
#                 pdf.output("Report.pdf")
  
    
        button3 = Button(app, text="Next", width=12, bg='#03A9F4', fg='#fff', command=getoutput)
        # button3.pack()
        button3.grid(row=23, column=0, padx=10, pady=10)
    button2 = Button(app, text="Check", width=12, bg='#03A9F4', fg='#fff', command=getcheck)
    # button2.pack()
    button2.grid(row=17, column=1, padx=10, pady=10)

# def getclear():
#     clear_span_entry.destroy()
#     width_Beam_entry.destroy()
#     width_support_entry.destroy()
#     clear_cover_entry.destroy()
#     overall_depth_entry.destroy()
#     dia_tension_reinforcement_entry.destroy()
#     dia_comp_reinforcement_entry.destroy()
#     dia_stirrups_entry.destroy()
#     grade_steel_entry.destroy()
#     grade_concrete_entry.destroy()
#     grade_stirrups_entry.destroy()
#     imposed_load_entry.destroy()
#     no_legs_stirrups_entry.destroy()
#
#     clear_span_entry = Entry(app, textvariable=clear_span_text, width=50);
#     clear_span_entry.grid(row=1, column=1);
#
#     width_Beam_entry = Entry(app, textvariable=width_Beam_text, width=50)
#     width_Beam_entry.grid(row=2, column=1)
#
#     width_support_entry = Entry(app, textvariable=width_support_text, width=50)
#     width_support_entry.grid(row=3, column=1)
#
#     clear_cover_entry = Entry(app, textvariable=clear_cover_text, width=50)
#     clear_cover_entry.grid(row=4, column=1)
#
#     overall_depth_entry = Entry(app, textvariable=overall_depth_text, width=50)
#     overall_depth_entry.grid(row=5, column=1)
#
#     dia_tension_reinforcement_entry = Entry(app, textvariable=dia_tension_reinforcement_text, width=50)
#     dia_tension_reinforcement_entry.grid(row=6, column=1)
#
#     dia_comp_reinforcement_entry = Entry(app, textvariable=dia_comp_reinforcement_text, width=50)
#     dia_comp_reinforcement_entry.grid(row=7, column=1)
#
#     dia_stirrups_entry = Entry(app, textvariable=dia_stirrups_text, width=50)
#     dia_stirrups_entry.grid(row=8, column=1)
#
#     grade_stirrups_entry = Entry(app, textvariable=grade_stirrups_text, width=50)
#     grade_stirrups_entry.grid(row=9, column=1)
#
#     grade_concrete_entry = Entry(app, textvariable=grade_concrete_text, width=50)
#     grade_concrete_entry.grid(row=10, column=1)
#
#     grade_steel_entry = Entry(app, textvariable=grade_steel_text, width=50)
#     grade_steel_entry.grid(row=11, column=1)
#
#     imposed_load_entry = Entry(app, textvariable=imposed_load_text, width=50)
#     imposed_load_entry.grid(row=12, column=1)
#
#     no_legs_stirrups_entry = Entry(app, textvariable=no_legs_stirrups_text, width=50)
#     no_legs_stirrups_entry.grid(row=13, column=1)
#

button1 = Button(app,text="Run",width=12,bg='#03A9F4',fg='#fff',command=getvals)
#button1.pack()
button1.grid(row=17,column=0,padx=10,pady=10)

# button4 = Button(app,text="Clear",width=12,bg='#03A9F4',fg='#fff',command=getclear)
# #button1.pack()
# button4.grid(row=15,column=0,padx=10,pady=10)

app.mainloop()


# In[ ]:





# In[ ]:




