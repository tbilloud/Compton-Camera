import sys
import pandas as pd
import uproot
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QPushButton, QWidget, QHBoxLayout,
    QTableWidget, QTableWidgetItem, QSplitter
)
from PyQt5.QtWebEngineWidgets import QWebEngineView
import plotly.express as px
import os  # For working with file paths
from PyQt5.QtCore import QUrl  # For handling URLs
from PyQt5.QtGui import QColor  # Import for setting background color


class MainWindow(QMainWindow):
    def __init__(self, df):
        super().__init__()
        self.df = df
        self.event_ids = df['EventID'].unique()
        self.current_event_index = 0
        self.table_visible = True
        self.removed_points = []  # Store removed points information

        self.initUI()

    def initUI(self):
        self.resize(1600, 600)
        self.setWindowTitle("3D Event Viewer with Optional Event Data Table")

        # Main layout and splitter
        main_widget = QWidget(self)
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)

        self.splitter = QSplitter(self)
        layout.addWidget(self.splitter)

        # Left side: Plotly-based 3D plot using QWebEngineView
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        plot_widget.setMinimumWidth(600)
        self.web_view = QWebEngineView()
        plot_layout.addWidget(self.web_view)
        self.splitter.addWidget(plot_widget)

        # Right side: Table for event data
        self.table_display = QTableWidget()
        self.table_display.setMinimumWidth(1000)
        self.splitter.addWidget(self.table_display)

        # Control buttons for navigation and table toggle
        button_layout = QHBoxLayout()
        self.prev_button = QPushButton("Previous Event")
        self.next_button = QPushButton("Next Event")
        self.toggle_table_button = QPushButton("Toggle Table")
        self.print_button = QPushButton("Print Event Data")
        button_layout.addWidget(self.prev_button)
        button_layout.addWidget(self.next_button)
        button_layout.addWidget(self.toggle_table_button)
        button_layout.addWidget(self.print_button)
        layout.addLayout(button_layout)

        # Connect buttons
        self.next_button.clicked.connect(self.next_event)
        self.prev_button.clicked.connect(self.prev_event)
        self.toggle_table_button.clicked.connect(self.toggle_table)
        self.print_button.clicked.connect(self.print_event_data)

        # Initial view updates
        self.plot_event()
        self.update_event_table()

    def create_plotly_figure(self, event_data):
        """Creates a Plotly 3D scatter plot for the event data with equalized axis ranges."""
        # Extract positions and other data
        x = event_data['Position_X']
        y = event_data['Position_Y']
        z = event_data['Position_Z']
        colors = event_data['ParticleName']

        # Change default range to have the same distance on all axes
        # Calculate ranges for each axis
        x_range = (x.min(), x.max())
        y_range = (y.min(), y.max())
        z_range = (z.min(), z.max())
        # Find the largest range among axes
        x_dist = x_range[1] - x_range[0]
        y_dist = y_range[1] - y_range[0]
        z_dist = z_range[1] - z_range[0]
        max_dist = max(x_dist, y_dist, z_dist)
        # Expand the ranges of smaller axes to match the largest range
        extra = 10 # add a little extra range to avoid points at the extremity to be cut by the viewer
        x_center = (x_range[0] + x_range[1]) / 2
        y_center = (y_range[0] + y_range[1]) / 2
        z_center = (z_range[0] + z_range[1]) / 2
        new_x_range = (x_center - max_dist / 2 - extra, x_center + max_dist / 2 + extra)
        new_y_range = (y_center - max_dist / 2 - extra, y_center + max_dist / 2 + extra)
        new_z_range = (z_center - max_dist / 2 - extra, z_center + max_dist / 2 + extra)

        # Create separate text annotations for gamma points
        text_labels = event_data['TrackID'].where(event_data['ParticleName'] == 'gamma', None)

        # Create 3D scatter plot
        fig = px.scatter_3d(
            event_data,
            x='Position_X',
            y='Position_Y',
            z='Position_Z',
            color=colors,
            title=f"Event ID: {self.event_ids[self.current_event_index]}",
            labels={'Position_X': 'X Position (um)', 'Position_Y': 'Y Position (um)', 'Position_Z': 'Z Position (um)'},
            text=text_labels  # This adds the `TrackID` as annotations for points with PDGEncoding == 22
        )

        fig.update_traces(marker=dict(size=5))
        fig.update_layout(
            legend_title_text="",
            scene=dict(
                xaxis=dict(title="Position_X", range=new_x_range),
                yaxis=dict(title="Position_Y", range=new_y_range),
                zaxis=dict(title="Position_Z", range=new_z_range)
            )
        )
        return fig

    def plot_event(self):
        """Render the 3D Plotly plot for the current EventID with PDGEncoding color mapping."""
        # Get current event data
        current_event_id = self.event_ids[self.current_event_index]
        event_data = self.df[self.df['EventID'] == current_event_id]

        event_data['Position_X'] = pd.to_numeric(event_data['Position_X'], errors='coerce')
        event_data['Position_Y'] = pd.to_numeric(event_data['Position_Y'], errors='coerce')
        event_data['Position_Z'] = pd.to_numeric(event_data['Position_Z'], errors='coerce')
        event_data['ParticleName'] = event_data['ParticleName'].astype(str)  # Ensure names are strings

        # Remove points that were deleted in the current event
        if self.removed_points:
            event_data = event_data[
                ~event_data.apply(
                    lambda row: (
                                    row['Position_X'],
                                    row['Position_Y'],
                                    row['Position_Z']
                                ) in self.removed_points,
                    axis=1,
                )
            ]

        # Create the Plotly figure
        fig = self.create_plotly_figure(event_data)

        # Get the figure as an HTML string
        plot_html = self.create_plotly_figure(event_data).to_html(include_plotlyjs='cdn')

        # Load the HTML string in QWebEngineView
        self.web_view.setHtml(plot_html)

    def update_event_table(self):
        """Displays rows corresponding to the current event in a table widget."""
        self.table_display.clear()

        # Get current event's data
        current_event_id = self.event_ids[self.current_event_index]
        event_data = self.df[self.df['EventID'] == current_event_id]

        # Drop unnecessary columns
        columns_to_drop = [
            'runID', 'time', 'stepLength', 'trackLength', 'posX', 'posY', 'posZ',
            'Position_X', 'Position_Y', 'Position_Z', 'nCrystalCompt', 'nCrystalRayl',
            'nCrystalConv', 'layerName', 'volumeID', 'trackLocalTime', 'sourcePosX',
            'sourcePosY', 'sourcePosZ', 'sourceEnergy', 'sourcePDG'
        ]
        cleaned_event_data = event_data.drop(columns=columns_to_drop, errors='ignore')

        # Update table structure
        self.table_display.setColumnCount(len(cleaned_event_data.columns))
        self.table_display.setRowCount(len(cleaned_event_data.index))
        self.table_display.setHorizontalHeaderLabels(cleaned_event_data.columns)

        # Populate table with data and highlight rows with PDGEncoding == 22
        for i, (index, row) in enumerate(cleaned_event_data.iterrows()):
            for j, value in enumerate(row):
                value = round(value, 3) if isinstance(value, float) else value
                item = QTableWidgetItem(str(value))

                # Highlight row if PDGEncoding == 22
                if row['ParticleName'] == 'gamma':
                    item.setBackground(QColor('yellow'))  # Set the background color to yellow

                self.table_display.setItem(i, j, item)

    def next_event(self):
        """Go to the next event."""
        if self.current_event_index < len(self.event_ids) - 1:
            self.current_event_index += 1
            self.plot_event()
            self.update_event_table()

    def prev_event(self):
        """Go to the previous event."""
        if self.current_event_index > 0:
            self.current_event_index -= 1
            self.plot_event()
            self.update_event_table()

    def toggle_table(self):
        """Show or hide the event data table."""
        if self.table_visible:
            self.splitter.widget(1).hide()
        else:
            self.splitter.widget(1).show()
        self.table_visible = not self.table_visible

    def print_event_data(self):
        """Print the DataFrame rows corresponding to the currently displayed event."""
        current_event_id = self.event_ids[self.current_event_index]
        event_data = self.df[self.df['EventID'] == current_event_id]
        print(event_data)


# Load ROOT file and initialize GUI
def main():
    # Read the ROOT file into a DataFrame
    file_path = "../output/"
    tree = uproot.open(file_path + "CC_Hits.root:Hits")
    df = tree.arrays(library="pd")

    # Unit conversions
    df['KineticEnergy'] *= 1000  # To keV
    df['TotalEnergyDeposit'] *= 1000  # To keV
    df['Position_X'] *= 1000  # To micrometers (um)
    df['Position_Y'] *= 1000  # To micrometers (um)
    df['Position_Z'] *= 1000  # To micrometers (um)

    # Launch PyQt5 app
    app = QApplication(sys.argv)
    main_window = MainWindow(df)
    main_window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
