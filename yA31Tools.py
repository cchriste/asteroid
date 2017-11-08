class HistExplorer:

    def __init__(self, Directory=None):
        if Directory==None:
            raise Exception("You must specify a subdirectory that contains all the yA31 timesteps in 300 cubed resolution.\nie- HistExplorer('/path/to/my/parent/directory/')")
        from vtk import vtkImageDataGeometryFilter, vtkPolyDataMapper, vtkXMLImageDataReader
        import numpy as np
        self.np = np
        self._Reader = vtkXMLImageDataReader
        self._load_files(Directory)

    def _load_files(self,Directory):
        import os
        self.Filenames = ['/'.join([Directory,i]) for i in os.listdir(Directory) if 'vti' in i]
        if len(self.Filenames)==0:
            raise Exception('Structured grid VTI files are exclusively supported in this version.')
        self.reader = self._Reader()
        self.reader.SetFileName(self.Filenames[0])
        self.reader.Update()
        self.my_step = 0;
        self.my_idx_max = self.reader.GetNumberOfPoints()

    # This function is quite slow. I am using a formula for structured grids that I have pulled from vtk source,
    #       instead of using this function. If a result seems fishy, consider using this function instead of 
    #       the direct formula for fact checking.
    # Source: GetLinearIndex at https://www.vtk.org/doc/nightly/html/vtkStructuredData_8h_source.html#l00250
    def _get_idx(self,my_3tuple):
        try:
            return(self.coord_map[my_3tuple][0])
        except AttributeError:
            try:
                if any(i<j for i,j in zip(self.getDimensions(),my_3tuple)):
                    raise Exception('Your array index is out of bounds')
                return(self.reader.GetOutput().ComputePointId(my_3tuple))
            except TypeError:
                raise Exception('Please enter a 3tuple within the dimensions given my HistExplorer().getDimensions()')

    def getSpacialCoord(self,my_3tuple):
        try:
            return(self.coord_map[my_3tuple][1])
        except AttributeError:
            return(self.reader.GetOutput().GetPoint(self._get_idx(my_3tuple)))
            
    def getTimesteps(self):
        return([idx for idx in range(len(self.Filenames))])

    def getStep(self):
        return(self.my_step)

    def nextStep(self):
        if self.my_step+1 > len(self.Filenames):
            raise Exception('You went out of bounds')
        self.reader.SetFileName(self.Filenames[self.my_step+1])
        self.reader.Update()
        self.my_step += 1;

    def lastStep(self):
        if self.my_step-1 < 0:
            raise Exception('You went out of bounds')
        self.reader.SetFileName(self.Filenames[self.my_step-1])
        self.reader.Update()
        self.my_step -= 1;

    def setStep(self,this_step):
        this_step = int(this_step)
        if this_step <= 0 or this_step > len(self.Filenames):
            raise Exception('You went out of bounds')
        self.reader.SetFileName(self.Filenames[this_step])
        self.reader.Update()
        self.my_step=this_step

    def getVars(self):
        try:
            self.my_vars
        except AttributeError:
            self.my_vars=[self.reader.GetOutput().GetPointData().GetArray(i).GetName() for i in range(self.reader.GetOutput().GetPointData().GetNumberOfArrays())]
        return(self.my_vars)

    def setVar(self,this_var):
        try:
            self.my_vars
        except AttributeError:
            self.getVars()
        try:
            self.this_var = (lambda name:[idx for idx,value in enumerate(self.my_vars) if value==name])(this_var)[0]
        except IndexError:
            raise Exception('You are attempting to set a variable that does not exist in these data. Explore available vars with HistExplorer().getVars()')

    def getActiveVar(self):
        try:
            return(self.my_vars[self.this_var])
        except AttributeError:
            raise Exception('You have not set an active variable. Set one with HistExplorer().setVar() or explore available vars with HistExplorer.().getVars()')

    def getDimensions(self):
        try:
            self.my_dims
        except AttributeError:
            self.my_dims = self.reader.GetOutput().GetDimensions()
        return(self.my_dims)

    def getValue(self,my_3tuple):
        if type(my_3tuple) is not tuple or len(my_3tuple) != len(self.getDimensions()):
            raise Exception('Specify a 3tuple within bounds given my HistExplorer().getDimensions()')
        try:
            _my_getter  = self.reader.GetOutput().GetPointData().GetArray(self.this_var).GetValue
            return(_my_getter(self._get_idx(my_3tuple)))
        except AttributeError:
            raise Exception('You have not set an active variable. Set one with HistExplorer().setVar() or explore available vars with HistExplorer.().getVars()')

    def getVarRange(self):
        try:
            self.this_var
        except AttributeError:
            self.getActiveVar()
        return(self.reader.GetOutput().GetPointData().GetArray(self.this_var).GetRange())
    def getAllValueDump(self, *args):
        if len(args) == 0:
            extent = range(self.my_idx_max)
            length = self.my_idx_max
        elif len(args) == 2:
            extent = args[0]
            length = args[1]
        else:
            raise Exception("HistExplorer().getAllValueDump() called incorrectly.")
        try:
            _my_getter  = self.reader.GetOutput().GetPointData().GetArray(self.this_var).GetValue
            extent_tmp = extent
            result_arr = self.np.empty(length)
            for idx, elem in enumerate(extent): result_arr[idx] = _my_getter(elem)
            return(result_arr)
        except AttributeError:
            raise Exception('You have not set an active variable. Set one with HistExplorer().setVar() or explore available vars with HistExplorer.().getVars()')

    def getExtentValues(self, *args):
        if len(args) != len(self.getDimensions()):
            raise Exception("You must specify as many range extents as you have dimensions.\nie- dim(300,300,300) -> HistExplorer()._get_tuples(range(200,300),range(100,150),range(300) is valid.\n    dim(300,300,300) -> HistExplorer()._get_tuples(range(100)) is not valid.")
        if len(args) != 3:
            raise Exception("This tool currently supports three dimensions")
        j_max = self.getDimensions()[1]
        i_max = self.getDimensions()[0]
        tuples = ((i,j,k) for i in args[0] for j in args[1] for k in args[2])
        my_extent_arr_idx = ((i,j,k) for i in range(len(args[0])) for j in range(len(args[1])) for k in range(len(args[2])))
        my_extent_idx = ((my_tuple[2]*j_max+my_tuple[1])*i_max+my_tuple[0] for my_tuple in tuples)
        result_arr = self.getAllValueDump(my_extent_idx,self.np.prod(self.getDimensions()))
        return_arr = self.np.zeros([len(i) for i in args])
        for this_idx,this_arr_idx in enumerate(my_extent_arr_idx):
            return_arr[this_arr_idx] = result_arr[this_idx]
        return(return_arr)

    # This is a wrapper/special case of getExtentValues()
    def getAllValues(self):
        try:
            return(self.getExtentValues(*[range(i) for i in self.getDimensions()]))
        except AttributeError:
            raise Exception('You have not set an active variable. Set one with HistExplorer().setVar() or explore available vars with HistExplorer.().getVars()')

    def getReader(self):
        return(self.reader)
