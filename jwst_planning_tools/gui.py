from jwst_planning_tools.util import get_v3PA_range
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class SimpleGUIApp:
    def __init__(self, root):
        self.root = root
        self.root.title("MRS Visibility Tool")

        # Initialize the GUI components
        self.create_widgets()

    def create_widgets(self):
        # Labels and input fields
        name_label = tk.Label(self.root, text="Name:")
        name_label.grid(row=0, column=0, padx=10, pady=5)

        self.name_entry = tk.Entry(self.root)
        self.name_entry.grid(row=0, column=1, padx=10, pady=5)
        self.name_entry.bind("<Return>", self.on_name_change)

        v3pa_label = tk.Label(self.root, text="V3PA:")
        v3pa_label.grid(row=1, column=0, padx=10, pady=5)

        self.v3pa_entry = tk.Entry(self.root)
        self.v3pa_entry.grid(row=1, column=1, padx=10, pady=5)
        self.v3pa_entry.bind("<Return>", self.plot_graph)
        # Message label
        self.message_label = tk.Label(self.root, text="Message will appear here.")
        self.message_label.grid(row=2, column=0, columnspan=2, padx=10, pady=10)


        plot_button = ttk.Button(self.root, text="Save Plot", command=self.save_plot())
        plot_button.grid(row=3, column=1, padx=10, pady=5)

    def on_name_change(self, event):
        out = get_v3PA_range(self.name_entry.get())
        if out is None:
            self.message = f"Name {self.name_entry.get()} could not be resolved by SIMBAD"
        else:
            start, end, ra, dec, v3pa_range = out
            self.v3pa_range = v3pa_range
            self.ra, self.dec = ra, dec
            self.message = f"{self.name_entry.get()}, has coordinates {self.ra}, {self.dec}. " \
                           f"\n Choose between the following V3PA ranges: {self.v3pa_range[0]} and {self.v3pa_range[1]}"
        self.update_message()

    def update_message(self):
        """Update the message label with the user's input."""
        name = self.name_entry.get()
        v3pa = self.v3pa_entry.get()
        self.message_label.config(text=self.message)

    def plot_graph(self):
        """Generate and display a simple plot based on the V3PA value."""
        fig, ax = plt.subplots()

        try:
            v3pa_value = int(self.v3pa_entry.get())
        except ValueError:
            self.message_label.config(text="Please enter a valid number for V3PA.")
            return

        # Dummy data for plotting
        x = [1, 2, 3, 4, 5]
        y = [i * v3pa_value for i in x]

        ax.plot(x, y, label="V3PA graph")
        ax.set_title("Sample Plot")
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        ax.legend()

        # Embed the plot in the Tkinter window
        canvas = FigureCanvasTkAgg(fig, master=self.root)
        canvas.draw()
        canvas.get_tk_widget().grid(row=5, column=0, columnspan=2)

# Run the application
if __name__ == "__main__":
    root = tk.Tk()
    app = SimpleGUIApp(root)
    root.mainloop()
