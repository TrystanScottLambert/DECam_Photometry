import tkinter

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import numpy as np

COUNTER = 0
POSSIBLE = []
ARTIFACT = []


def start_gui(i_bands, z_bands, n_bands):

    root = tkinter.Tk()
    root.wm_title("Finding Artifacts")

    fig = Figure(figsize=(5, 4), dpi=100)
    t = np.arange(0, 3, .01)
    ax_i = fig.add_subplot(131)
    i_band_display = ax_i.imshow(i_bands[COUNTER])
    ax_i.set_xlabel('i-band')
    ax_i.set_xticks([])
    ax_i.set_yticks([])

    ax_z = fig.add_subplot(132)
    z_band_display = ax_z.imshow(z_bands[COUNTER])
    ax_z.set_xlabel('z-band')
    ax_z.set_xticks([])
    ax_z.set_yticks([])

    ax_n = fig.add_subplot(133)
    n_band_display = ax_n.imshow(n_bands[COUNTER])
    ax_n.set_xlabel('n-band')
    ax_n.set_xticks([])
    ax_n.set_yticks([])

    canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
    canvas.draw()

    # pack_toolbar=False will make it easier to use a layout manager later on.
    toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
    toolbar.update()

    canvas.mpl_connect(
        "key_press_event", lambda event: print(f"you pressed {event.key}"))
    canvas.mpl_connect("key_press_event", key_press_handler)

    button_quit = tkinter.Button(master=root, text="Quit", command=root.destroy)


    def update_frequency():
        global COUNTER
        COUNTER += 1

        i_band_display.set_data(i_bands[COUNTER])
        z_band_display.set_data(z_bands[COUNTER])
        n_band_display.set_data(n_bands[COUNTER])

        canvas.draw()


    def update_possible():
        global POSSIBLE
        global COUNTER
        POSSIBLE.append(COUNTER)
        update_frequency()

    def update_artifact():
        global ARTIFACT
        global COUNTER
        ARTIFACT.append(COUNTER)
        update_frequency()


    button_no = tkinter.Button(root, command=update_artifact, text='REJECT',bg='red')
    button_yes = tkinter.Button(root, command=update_possible, text='MAYBE',bg='green')
    # Packing order is important. Widgets are processed sequentially and if there
    # is no space left, because the window is too small, they are not displayed.
    # The canvas is rather flexible in its size, so we pack it last which makes
    # sure the UI controls are displayed as long as possible.
    button_quit.pack(side=tkinter.BOTTOM)
    button_yes.pack(side=tkinter.TOP)
    button_no.pack(side=tkinter.TOP)
    toolbar.pack(side=tkinter.BOTTOM, fill=tkinter.X)
    canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

    tkinter.mainloop()
    return ARTIFACT, POSSIBLE
