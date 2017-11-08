from visuspy import *

class VtuToIdx:
  
  def testIdx(self):
    self.filename="temp/tutorial_1.idx"
    self.WriteIdx()
    self.ReadIdx()
    self.MergeIdx()
  
  # WriteIdx
  def WriteIdx(self): 
    
    dataset_box=NdBox(NdPoint(0,0,0),NdPoint.one(16,16,16))
    
    idxfile=IdxFile();
    idxfile.box=NdBox(dataset_box)
    idxfile.fields.push_back(Field("myfield",DType.parseFromString("uint32")))

    bSaved=idxfile.save(self.filename)
    self.assertTrue(bSaved)
    
    dataset=Dataset.loadDataset(self.filename)
    self.assertIsNotNone(dataset)
    access=dataset.get().createAccess()
    
    sampleid=0
    
    for Z in range(0,16):
      slice_box=dataset.get().getBox().getZSlab(Z,Z+1)
      
      query=QueryPtr(Query(dataset.get(),ord('w')))
      query.get().position=Position(slice_box)
      
      self.assertTrue(dataset.get().beginQuery(query))
      self.assertEqual(query.get().nsamples.innerProduct(),16*16)
      
      buffer=Array(query.get().nsamples,query.get().field.dtype)
      query.get().buffer=buffer
      
      fill=convertToNumPyArray(buffer)
      for Y in range(16):
        for X in range(16):
          fill[Y,X]=sampleid
          sampleid+=1

      self.assertTrue(dataset.get().executeQuery(access,query))

  # ReadIdx
  def ReadIdx(self): 
    
    dataset=Dataset_loadDataset(self.filename)
    self.assertIsNotNone(dataset)
    box=dataset.get().getBox()
    field=dataset.get().getDefaultField()
    access=dataset.get().createAccess()
    
    sampleid=0
    for Z in range(0,16):
      slice_box=box.getZSlab(Z,Z+1)
      
      query=QueryPtr(Query(dataset.get(),ord('r')))
      query.get().position=Position(slice_box)
      
      self.assertTrue(dataset.get().beginQuery(query))
      self.assertEqual(query.get().nsamples.innerProduct(),16*16)
      self.assertTrue(dataset.get().executeQuery(access,query))
      
      check=convertToNumPyArray(query.get().buffer)
      for Y in range(16):
        for X in range(16):
          self.assertEqual(check[Y,X],sampleid)
          sampleid+=1

  def MergeIdx(self): 
    
    dataset=Dataset_loadDataset(self.filename)
    self.assertIsNotNone(dataset)
    
    box=dataset.get().getBox()
    access=dataset.get().createAccess()
    field=dataset.get().getDefaultField()
    MaxH=dataset.get().getBitmask().getMaxResolution()
    self.assertEqual(MaxH,12) #in the bitmask_pattern "V012012012012" the very last bit of the bitmask is at position MaxH=12 
    
    #I want to read data from first slice Z=0
    slice_box=box.getZSlab(0,1);
    
    #create and read data from VisusFIle up to resolution FinalH=8 (<MaxH)
    query=QueryPtr(Query(dataset.get(),ord('r')))
    query.get().position=Position(slice_box)
    query.get().end_resolutions.push_back(8)
    query.get().end_resolutions.push_back(12)
    query.get().merge_mode=Query.InsertSamples
    
    # end_resolution=8
    
    self.assertTrue(dataset.get().beginQuery(query))
    self.assertTrue(query.get().nsamples.innerProduct()>0)
    self.assertTrue(dataset.get().executeQuery(access,query))
    self.assertEqual(query.get().cur_resolution,8)
    
    # end_resolution=12
    self.assertTrue(dataset.get().nextQuery(query))
    self.assertEqual(query.get().nsamples.innerProduct(),16*16)
    self.assertTrue(dataset.get().executeQuery(access,query))
    self.assertEqual(query.get().cur_resolution,12)
    
    #verify the data is correct
    check=convertToNumPyArray(query.get().buffer)
    sampleid=0
    for Y in range(0,16):
      for X in range(0,16):
        self.assertEqual(check[Y,X],sampleid)
        sampleid+=1 
        
    # finished
    self.assertFalse(dataset.get().nextQuery(query)) 


class VtuUtils:
    
    def __init__(self, base_filename="/usr/sci/cedmav/data/asteroid/pv_insitu_39966/pv_insitu_39966_0_", nparts=1024
                 ,idx_filename='/usr/sci/cedmav/data/asteroid/converted/idx/asteroid_full.idx'):
        from vtk import vtkXMLUnstructuredGridReader
        self.vtkXMLUnstructuredGridReader = vtkXMLUnstructuredGridReader
        self.base_filename=base_filename
        self.idx_filename=idx_filename
        self.nparts=nparts
        self.bounds_min=[1e12,1e12,1e12]
        self.bounds_max=[-1e12,-1e12,-1e12]
        self.amr_counts=[0,0,0,0]  #assumes four amr levels
        SetCommandLine()
        IdxModule.attach()

        self.index_shift=[2300000.0, 500000.0, 1200000.0]   #note: hardcoded, but this is what getbounds returned
        self.voxel_delta=2500

    def getFieldNames(self):
        rdr=self.vtkXMLUnstructuredGridReader()
        filename=self.base_filename+str(0)+".vtu"
        rdr.SetFileName(filename)
        rdr.Update()
        celldata=rdr.GetOutput().GetCellData()

        self.fields=[]
        for i in range(celldata.GetNumberOfArrays()):
            self.fields.append(celldata.GetArrayName(i))

        print "**** Fields:"
        print self.fields

    def getBoundsAndNumCells(self):
        rdr=self.vtkXMLUnstructuredGridReader()
        print "Getting bounds and amr cell counts..."

        update_frequency=10;
        update_output=0;
        for i in range(0,self.nparts):
            filename=self.base_filename+str(i)+".vtu"
            rdr.SetFileName(filename)
            rdr.Update()
            pdata=rdr.GetOutput().GetPoints().GetData()
            #print "checking "+str(pdata.GetNumberOfTuples())+" points in "+filename
            for p in range(pdata.GetNumberOfTuples()):
                X=pdata.GetTuple3(p)
                self.bounds_min[0]=min(self.bounds_min[0],X[0])
                self.bounds_max[0]=max(self.bounds_max[0],X[0])
                self.bounds_min[1]=min(self.bounds_min[1],X[1])
                self.bounds_max[1]=max(self.bounds_max[1],X[1])
                self.bounds_min[2]=min(self.bounds_min[2],X[2])
                self.bounds_max[2]=max(self.bounds_max[2],X[2])

            arr=rdr.GetOutput().GetCellData().GetArray('grd')
            #print "checking "+str(arr.GetNumberOfValues())+" cells in "+filename
            for k in range(arr.GetNumberOfValues()):
                val=arr.GetValue(k)
                self.amr_counts[int(val)-1]+=1

            if update_output >= update_frequency:
                print "Current range (",str(i),"/",self.nparts,"blocks checked):",self.bounds_min,"to",self.bounds_max
                print "Current counts (",str(i),"/",self.nparts,"blocks checked):",self.amr_counts
                update_output=0
            update_output+=1

        print "**** Final Range:"
        print self.bounds_min,"to",self.bounds_max
        print "**** Final Counts:"
        print self.amr_counts

    def createIdx(self):
        self.getFieldNames()

        # assumes fields are float32
        fields=self.fields
        fields.remove('grd')  #this is just the amr level

        #logic_to_physic=NdBox(NdPoint(self.bounds_min[0],self.bounds_min[1],self.bounds_min[2]),NdPoint.one(self.bounds_max[0],self.bounds_max[1],self.bounds_max[2]))
        dataset_box=NdBox(NdPoint(0,0,0),NdPoint.one(1840,1120,960))
        idxfile=IdxFile();
        idxfile.box=NdBox(dataset_box)
        for field in fields:
            idxfile.fields.push_back(Field(field,DType.parseFromString("float32")))
        idxfile.timesteps=DatasetTimesteps(0,477)
        idxfile.time_template="time%05d/"

        bSaved=idxfile.save(self.idx_filename)
        if bSaved:
            print "idx file created successfully!"
        else:
            print "error creating idx file"
            

    def testOpenIdx(self):
        dataset=Dataset.loadDataset(self.idx_filename)
        #self.assertIsNotNone(dataset)
        access=dataset.get().createAccess()
    
    def realToIndex(self,pt):
        return [int((self.index_shift[0]+pt[0])/self.voxel_delta),
                int((self.index_shift[1]+pt[1])/self.voxel_delta),
                int((self.index_shift[2]+pt[2])/self.voxel_delta)]


    def writeIdx(self,fieldname,start_part=0):
        from vtk import vtkIdList
        
        dataset=Dataset.loadDataset(self.idx_filename)
        max_resolution=dataset.get().getMaxResolution()
        print "dataset max resolution is",max_resolution
        access=dataset.get().createAccess()
        field=dataset.get().getFieldByName(fieldname)

        #want to reuse a query (and its buffers), but no way to reset status once it's been executed
        #query=QueryPtr(Query(dataset.get(),ord('w')))
        #query.get().field=field

        ##test query
        # import pdb
        # idxpt=[0,0,0]
        # query_box=NdBox(NdPoint(idxpt[0],idxpt[1],idxpt[2]),NdPoint.one(idxpt[0]+1,idxpt[1]+1,idxpt[2]+1))
        # query_box_from_dataset=dataset.get().getBox().getZSlab(0,1)
        # print "query box:",query_box.toString()
        # query.get().position=Position(query_box)
        # ret=dataset.get().beginQuery(query)
        # if not ret:
        #     print "begin query failed"
        # print "query size:",query.get().nsamples.innerProduct()
        # print "query nsamples:",query.get().nsamples.toString()
        # pdb.set_trace()

        # arr=Array(query.get().nsamples,query.get().field.dtype)
        # print "array dims:",arr.dims.toString()
        # query.get().buffer=arr
        # fill=convertToNumPyArray(arr).reshape((1))
        # fill[0]=42.0

        # ret=dataset.get().executeQuery(access,query)
        # # if not ret:
        # #     print "write query failed"
        # if ret:
        #     print "write query succeeded"
        # return
        ##test query

        import pdb
        rdr=self.vtkXMLUnstructuredGridReader()
        for i in range(start_part,self.nparts):
            filename=self.base_filename+str(i)+".vtu"
            rdr.SetFileName(filename)
            rdr.Update()
            output=rdr.GetOutput()
            pdata=output.GetPoints().GetData()
            celldata=output.GetCellData()
            amr_level=celldata.GetArray('grd')
            field_data=celldata.GetArray(fieldname)
            pts=vtkIdList()
            update_frequency=10000;
            update_output=0;

            for cell_id in range(0,output.GetNumberOfCells()):
                output.GetCellPoints(cell_id,pts)
                #todo: assert that first point is always the "smallest" (closest to bounds_min)
                #for idx in range(pts.GetNumberOfIds()):
                    #todo

                query=QueryPtr(Query(dataset.get(),ord('w')))
                query.get().field=field
                set_resolution=max_resolution - 2*(4-int(amr_level.GetValue(cell_id)))
                #if set_resolution < max_resolution:
                #    print "setting query max resolution to",set_resolution
                #query.get().max_resolution = set_resolution   #no! the query has the same max_res as the dataset, but needs to stop before that using end_resolutions
                query.get().end_resolutions.push_back(set_resolution)
                #print "max_resolution:",query.get().max_resolution

                pt=pdata.GetTuple3(pts.GetId(0))  #assuming the first is the "root" of the voxel
                idxpt=self.realToIndex(pt)
                query_box=NdBox(NdPoint(idxpt[0],idxpt[1],idxpt[2]),NdPoint.one(idxpt[0]+1,idxpt[1]+1,idxpt[2]+1))
                #print "query box:",query_box.toString()
                query.get().position=Position(query_box)
       
                ret=dataset.get().beginQuery(query)
                if not ret:
                    print "begin query failed"
                    pdb.set_trace()
                #print "query size:",query.get().nsamples.innerProduct()
                #print "query nsamples:",query.get().nsamples.toString()
                #pdb.set_trace()
      
                arr=Array(query.get().nsamples,query.get().field.dtype)
                #print "array dims:",arr.dims.toString()
                query.get().buffer=arr
                fill=convertToNumPyArray(arr).reshape((1))
                fill[0]=field_data.GetValue(cell_id)

                ret=dataset.get().executeQuery(access,query)
                if not ret:
                    print "write query failed"
                    pdb.set_trace()
                #if ret:
                #    print "write query succeeded"

                if update_output >= update_frequency:
                    print "Written ",cell_id,"/",output.GetNumberOfCells(),"cells of this file"
                    update_output=0
                update_output+=1
    
            print "Finished with file",str(i),"/",self.nparts
        

# ////////////////////////////////////////////////////////
if __name__ == '__main__':
  SetCommandLine()
  IdxModule.attach()

  # do something...
  print "doing nothing..."

  IdxModule.detach()


