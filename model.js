<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Py3Dmol Example</title>
    <!-- Include Py3Dmol library -->
    <script src="https://3dmol.org/build/3Dmol.js"></script>
</head>
<body>
    <div id="pdbViewer" style="width: 800px; height: 600px;"></div>

    <script>
        // Create a viewer instance
        var viewer = $3Dmol.createViewer('pdbViewer', {backgroundColor: 'white'});

        // Load the PDB data from a URL or file
        var pdbData = `{{'1crn.pdb'}}`;

        // Add the PDB data to the viewer
        viewer.addModel(pdbData, 'pdb');

        // Set the style for the structure
        viewer.setStyle({}, {cartoon: {color: 'spectrum'}});

        // Zoom to fit the structure
        viewer.zoomTo();

        // Render the viewer
        viewer.render();
    </script>
</body>
</html>