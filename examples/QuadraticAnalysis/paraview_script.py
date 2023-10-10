from paraview.simple import *


def create_contour_from_vts(file_path, contour_value, density_name):
    # Read the VTS file
    reader = XMLStructuredGridReader(FileName=[file_path])

    # Update to make sure the file is read (important for distributed setups)
    UpdatePipeline()

    # Create contour
    contour = Contour(Input=reader)
    contour.ContourBy = [density_name]
    contour.Isosurfaces = [contour_value]

    # Show the contour
    Show(contour)

    # Update to generate the contour
    Render()


# Set the file path to your .vts file
file_path = "/path/to/your/file.vts"

# Create contour for ground-state density with contour_value=0.5
create_contour_from_vts(file_path, 0.5, "ground_density")

# Create contour for response density with contour_value=0.5
create_contour_from_vts(file_path, 0.5, "response_density")
