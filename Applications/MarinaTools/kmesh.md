## *kmesh* Usage
* -o
	* Specify a filename for output; used with other options
	* *ex)* -o filename.nrrd
* -scalarName
	* scalar name [string]
* -outputScalarName
	* scalar name for output [string]
* -sigma
	* sigma value [double]
* -threshold
	* Threshold value [double]
* -iter
	* number of iterations [int]
* -attrDim
	* The number of components of attribute
	* *ex)* -attrDim 3 (vector)
* -exportScalars
	* Export scalar values to a text file
	* *ex)* -exportScalars [in-mesh] [scalar.txt]
* -importScalars
	* Add scalar values to a mesh [in-mesh] [scalar.txt] [out-mesh]
* -smoothScalars
	* Gaussian smoothing of scalar values of a mesh. [in-mesh] [out-mesh]
* -copyScalars
	* Copy a scalar array of the input model to the output model
	* *ex)* -copyScalars input-model1 input-model2 output-model -scalarName name
* -averageScalars
	* Compute the average of scalars across given inputs
	* *ex)* -averageScalars -o output-vtk input1-vtk input2-vtk ... 
* -sampleImage
	* Sample pixels for each point of a given model. Currently, only supported image type is a scalar
	* *ex)* -sampleImage image.nrrd model.vtp output.vtp -outputScalarName scalarName
* -voronoiImage
	* Compute the voronoi image from a given data set. A reference image should be given.
	* *ex)* -voronoiImage ref-image.nrrd input-dataset output-image.nrrd -scalarName voxelLabel
* -scanConversion
	* Compute a binary image from a surface model
	* *ex)* -scanConversion input-surface input-image.nrrd output-image.nrrd
* -appendData
	* Append input meshes into a single data [output-mesh]
* -computeCurvature
	* Compute curvature values for each point
	* *ex)* -computeCurvature input-vtk output-vtk
* -vti
	* Convert an ITK image to VTI format (VTKImageData)
	* *ex)* -vti imageFile outputFile [-attrDim 3] [-maskImage mask]
* -vtu
	* Convert an ITK image to VTU format (vtkUnstructuredGrid). This is useful when masking is needed.
	* *ex)* -vtu imageFile outputFile -maskImage maskImage
* -maskImage
	* A mask image for the use of -vtu
	* *ex)* -maskImage mask.nrrd
* -traceStream
	* Trace a stream line from a given point set
	* *ex)* -traceStream input-vtu-field input-vtk output-lines output-points
* -traceDirection
	* Choose the direction of stream tracing (both, forward, backward)
	* *ex)* -traceStream ... -traceDirection (both|forward|backward)
* -zrotate
	* Rotate all the points along the z-axis. Change the sign of x and y coordinate.
	* *ex)* -traceStream ... -zrotate
* -traceClipping
	* Clip stream lines to fit with an object
	* *ex)* -traceClipping stream_lines.vtp stream_object.vtp stream_lines_output.vtp
* -traceScalarCombine
	* Combine scalar values from a seed object to a stream line object. The stream line object must have PointIds for association. -zrotate option will produce the rotated output.
	* *ex)* -traceScalarCombine stream_seed.vtp stream_lines.vtp stream_lines_output.vtp -scalarName scalarToBeCopied
* -filterStream
	* Filter out stream lines which are lower than a given threshold
	* *ex)* -filterStream stream-line-input stream-seed-input stream-line-output -scalarName scalar -threshold xx
* -thresholdMin
	* Give a minimum threshold value for -filterStream
	* *ex)* -threshold 10 (select a cell whose attriubte is greater than 10)
* -thresholdMax
	* Give a maximum threshold value for -filterStream
	* *ex)* -threshold 10 (select a cell whose attriubte is lower than 10)
* -fitting
	* Fit a model into a binary image
	* *ex)* -fitting input-model binary-image output-model
* -ellipse
	* Create an ellipse with parameters []
	* *ex)* -ellipse 101 101 101 51 51 51 20 20 20 -o ellipse.nrrd
* -h
	* print help message
