import sys, os, subprocess, fileinput, math, tempfile, datetime

class JobConfig():
    '''Class for storing configuration for each dataset'''
    def __init__(self, dataName = 'TEST', inDir = '/uscms/home/naodell/nobackup/TEST', nJobs = 1, arguments = 'TEST 0 muon'):
        self._dataName  = dataName
        self._inDir     = inDir
        self._nJobs     = nJobs
        self._args      = arguments

class BatchMaster():
    '''A tool for submitting batch jobs'''
    def __init__(self, configList, shortQueue = False, stageDir = 'outputBatch', executable = 'execBatch.csh', selection = 'test'):
        self._selection  = selection
        self._current    = os.path.abspath('.')
        self._stageDir   = stageDir
        self._configList = configList
        self._executable = executable
        self._shortQueue = shortQueue
    

    def get_current_time(self):
        ''' 
        Returns a string of the current time with
        the format  
        '''

        now = datetime.datetime.now()
        currentTime = '{0:02d}{1:02d}{2:02d}_{3:02d}{4:02d}{5:02d}'.format(now.year, now.month, now.day, now.hour, now.minute, now.second)
        return currentTime


    def make_directory(self, filePath, clear = True):
        '''
        Create save path in case it doesn't already exist
        '''

        if not os.path.exists(filePath):
            os.system('mkdir -p '+filePath)

        if clear and len(os.listdir(filePath)) != 0:
            os.system('rm '+filePath+'/*')


    def split_jobs(self, directory, nJobs):
        '''
        Split jobs by dataset
        '''

        fileList = os.listdir(directory)
        nFiles = len(fileList)
        
        # Split files to requested number.  Cannot exceed
        # the number of files being run over.
        if nJobs > nFiles:
            nJobs = nFiles

        nFilesPerJob = int(math.ceil(float(nFiles)/float(nJobs)))
        fileSplit = [fileList[i:i+nFilesPerJob] for i in range(0, len(fileList), nFilesPerJob)]


        return fileSplit


    def make_batch_lpc(self, cfg, count, sourceFiles):
        '''
        Prepares for submission to lpc.  Does the following:

        1. Generates input_files.txt with files to run over
        2. Write batch configuration file
        4. ????
        5. Poor house :(
        '''

        ## make file with list of inputs ntuples for the analyzer
        input = open('input_{1}_{2}.txt'.format(self._stageDir, cfg._dataName, str(count+1)), 'w')

        path = cfg._inDir
        if cfg._inDir[:5] == '/pnfs':
            path = 'dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr' + cfg._inDir[5:]

        for i,source in enumerate(sourceFiles):
            input.write(path+'/'+source+'\n')
        input.close()

        ## Writing the batch config file
        batch_tmp = open('.batch_tmp_{1}_{2}'.format(self._stageDir, cfg._dataName, count+1), 'w')
        batch_tmp.write('Universe              = vanilla\n')

        ## TESTING CONDORG ## --> Doesn't seem to work :(
        #batch_tmp.write('Universe              = globus\n')
        #batch_tmp.write('globusscheduler       = cmsosgce.fnal.gov/jobmanager-condor\n')

        batch_tmp.write('Arguments             = {0} {1} {2}\n'.format(cfg._dataName, count+1, cfg._args))
        batch_tmp.write('Executable            = {0}\n'.format(self._executable))
        #batch_tmp.write('Executable            = exe\n')
        batch_tmp.write('Should_Transfer_Files = YES\n')
        batch_tmp.write('WhenToTransferOutput  = ON_EXIT\n')
        batch_tmp.write('Transfer_Input_Files  = source.tar.gz, input_{1}_{2}.txt\n'.format(self._stageDir, cfg._dataName,  str(count+1)))
        batch_tmp.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
        batch_tmp.write('Notification          = Never\n')
        if self._shortQueue:
            batch_tmp.write('+LENGTH               = "SHORT"\n')
        batch_tmp.write('Output                = reports/'+cfg._dataName+'_'+str(count)+'_$(Cluster)_$(Process).stdout\n')
        batch_tmp.write('Error                 = reports/'+cfg._dataName+'_'+str(count)+'_$(Cluster)_$(Process).stderr\n')
        batch_tmp.write('Log                   = reports/'+cfg._dataName+'_'+str(count)+'_$(Cluster)_$(Process).log   \n')
        batch_tmp.write('notify_user           = brian.lee.pollack@cern.ch\n')
        batch_tmp.write('Queue\n')
        batch_tmp.close()

    def make_batch_lpc_alt(self, cfg, sources):
        '''
        Prepares for submission to lpc.  Does the following:

        1. Generates input_files.txt with files to run over
        2. Write batch configuration file
        4. ????
        5. Poor house :(
        '''

        ## Writing the batch config file
        batch_tmp = open('.batch_tmp_{0}'.format(cfg._dataName,), 'w')
        batch_tmp.write('Universe              = vanilla\n')
        batch_tmp.write('Should_Transfer_Files = YES\n')
        batch_tmp.write('WhenToTransferOutput  = ON_EXIT\n')
        batch_tmp.write('Notification          = Never\n')
        batch_tmp.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
        batch_tmp.write('notify_user           = brian.lee.pollack@cern.ch\n')
        if self._shortQueue:
            batch_tmp.write('+LENGTH               = "SHORT"\n')
        batch_tmp.write('\n')

        for i, source in enumerate(sources):
            ## make file with list of inputs ntuples for the analyzer
            input = open('input_{1}_{2}.txt'.format(self._stageDir, cfg._dataName, str(i+1)), 'w')

            path = cfg._inDir
            if cfg._inDir[:5] == '/pnfs':
                path = 'dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr' + cfg._inDir[5:]

            for file in source:
                input.write(path+'/'+file+'\n')
            input.close()

            batch_tmp.write('Arguments             = {0} {1} {2}\n'.format(cfg._dataName, i+1, cfg._args))
            batch_tmp.write('Executable            = {0}\n'.format(self._executable))
            batch_tmp.write('Transfer_Input_Files  = source.tar.gz, input_{0}_{1}.txt\n'.format(cfg._dataName, i+1))
            batch_tmp.write('Output                = reports/{0}_{1}_$(Cluster)_$(Process).stdout\n'.format(cfg._dataName, i+1))
            batch_tmp.write('Error                 = reports/{0}_{1}_$(Cluster)_$(Process).stderr\n'.format(cfg._dataName, i+1))
            batch_tmp.write('Log                   = reports/{0}_{1}_$(Cluster)_$(Process).log   \n'.format(cfg._dataName, i+1))
            batch_tmp.write('Queue\n\n')

        batch_tmp.close()
        

    def submit_to_batch(self, bSystem = 'lpc'):
        '''
        Submits batch jobs to batch.  Currently only works
        for lpc batch system, but should be updated for more 
        general use
        '''

        self._stageDir  = self._stageDir + '/' + self.get_current_time()
        self.make_directory(self._stageDir, clear=False)

        ## Creating tarball of current workspace
        os.system('tar czf {0}/source.tar.gz ../../analysis 2> /dev/null'.format(self._stageDir))

        subprocess.call('cp {0} {1}'.format(self._executable, self._stageDir), shell=True)
        os.chdir(self._stageDir)
        self.make_directory('reports', clear=False)
        
        if bSystem is 'lpc':

            for cfg in self._configList:
                sourceFiles = self.split_jobs(cfg._inDir, cfg._nJobs)

                #for i, source in enumerate(sourceFiles):
                #    self.make_batch_lpc(cfg, i, source)
                #    subprocess.call('condor_submit .batch_tmp_{1}_{2}'.format(self._stageDir, cfg._dataName, i+1), shell=True)

                self.make_batch_lpc_alt(cfg, sourceFiles)
                subprocess.call('condor_submit .batch_tmp_{0}'.format(cfg._dataName), shell=True)


