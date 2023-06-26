import tkinter

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
from astropy.visualization import ZScaleInterval

COUNTER = 0
POSSIBLE = []
ARTIFACT = []


def get_limits(array_2d: np.ndarray) -> tuple[float, float]:
    """
    Returns the Zscale limits for the given data.
    """
    zscale = ZScaleInterval()
    v_min, v_max = zscale.get_limits(array_2d)
    return v_min, v_max

def start_gui(i_bands, z_bands, n_bands):

    root = tkinter.Tk()
    root.wm_title("Finding Artifacts")

    fig = Figure(figsize=(5, 4), dpi=100)
    ax_i = fig.add_subplot(131)
    i_data = i_bands[COUNTER]
    i_limits = get_limits(i_data)

    i_band_display = ax_i.imshow(i_data, vmin = i_limits[0], vmax = i_limits[1])
    ax_i.set_xlabel('i-band')
    ax_i.set_xticks([])
    ax_i.set_yticks([])

    ax_z = fig.add_subplot(132)
    z_data = z_bands[COUNTER]
    z_limits = get_limits(z_data)
    z_band_display = ax_z.imshow(z_data, vmin=z_limits[0], vmax=z_limits[1])
    ax_z.set_xlabel('z-band')
    ax_z.set_xticks([])
    ax_z.set_yticks([])

    ax_n = fig.add_subplot(133)
    n_data = n_bands[COUNTER]
    n_limits = get_limits(n_data)
    n_band_display = ax_n.imshow(n_data, vmin=n_limits[0], vmax=n_limits[1])
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
        print(COUNTER, len(i_bands))
        POSSIBLE.append(COUNTER)
        if COUNTER == len(i_bands) - 1:
            root.destroy()
        else:
            update_frequency()

    def update_artifact():
        global ARTIFACT
        global COUNTER
        print(COUNTER, len(i_bands))
        ARTIFACT.append(COUNTER)
        if COUNTER == len(i_bands) - 1:
            root.destroy()
        else:
            update_frequency()


    button_no = tkinter.Button(root, command=update_artifact, text='REJECT',bg='red')
    button_yes = tkinter.Button(root, command=update_possible, text='MAYBE',bg='green')

    button_quit.pack(side=tkinter.BOTTOM)
    button_yes.pack(side=tkinter.TOP)
    button_no.pack(side=tkinter.TOP)
    toolbar.pack(side=tkinter.BOTTOM, fill=tkinter.X)
    canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

    tkinter.mainloop()
    return ARTIFACT, POSSIBLE
