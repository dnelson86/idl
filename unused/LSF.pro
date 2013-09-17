; LSF.pro
; automation and queue job related (pre-SLURM)
; dnelson jun.2013

; makeArepoFoFBsub(): write bsub file to invoke FoF/Subfind group finder post-processing

pro makeArepoFoFBsub

  ; config
  res = 512
  run = 'gadget'
  ;f = '10'
  
  sP = simParams(res=res,run=run)

  snapZ = redshiftToSnapNum([6.0,5.0,4.0,3.0,2.0,1.0,0.0],sP=sP)
  print,snapZ

  snapRange = [217,217,1]
  ;snaps = [257]

  ; job config
  spawnJobs = 1 ; execute bsub?
  nProcs    = 64 ; needed nodes: 128^3=0.5 (n4tile4 PartAllocFactor=2)
                 ;               256^3=4 (n32tile4 PartAllocFactor=2.5 MaxMemSize=7900)
                 ;               512^3=16 (n64tile4 PartAllocFactor=1.5 though n32tile4 ok until z=2)
  ptile     = 4 ; span[ptile=X]
  cmdCode   = 3 ; fof/substructure post process
  
  excludeHosts = ['hero1015','hero0103','hero2205','hero0410','hero2413','hero1004'] ;leave empty otherwise
 
  paramFile = "param_fof.txt"

  for snap=snapRange[0],snapRange[1],snapRange[2] do begin
  ;foreach snap,snaps do begin
  
    ; write bjob file
    jobFileName = sP.plotPath + 'job_fof.bsub'
    
    ; check before overriding
    if file_test(jobFileName) then message,'Error: job_fof.bsub already exists'
    
    print,'['+str(snap)+'] Writing: '+jobFilename
    
    ; host exclusion string
    selectStr = ''
    if n_elements(excludeHosts) gt 0 then begin
      selectStr = 'select['
      for i=0,n_elements(excludeHosts)-1 do selectStr += 'hname!='+excludeHosts[i]+' && '
      selectStr = strmid(selectStr,0,strlen(selectStr)-4) + '] '
    endif
      
    openW, lun, jobFilename, /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#BSUB -q nancy'
    printf,lun,'#BSUB -J fof_' + str(snap) + ''
    printf,lun,'#BSUB -n ' + str(nProcs)
    printf,lun,'#BSUB -R "' + selectStr + 'span[ptile=' + str(ptile) + ']"' ;rusage[mem=31000]
    printf,lun,'#BSUB -x'
    printf,lun,'#BSUB -o run_fof_'+str(snap)+'.out'
    printf,lun,'#BSUB -g /dnelson/fof' ; 4 concurrent jobs limit automatically
    printf,lun,'#BSUB -cwd ' + sP.arepoPath
    printf,lun,'#BSUB -a openmpi'
    printf,lun,''
      
    ; write projection commands
    strArray = ['mpirun.lsf ./Arepo_fof '+paramFile,$
                str(cmdCode),$
                str(snap),$
                '>> run_fof_'+str(snap)+'.txt']
    printf,lun,strjoin(strArray,' ')
      
    ; close
    close,lun
    free_lun,lun
      
    ; add to queue if requested
    if (spawnJobs) then begin
      spawn, 'bsub < ' + sP.plotPath + 'job_fof.bsub', result
      print,'  '+result
      
      ; file cleanup
      wait,0.5
      spawn, 'rm ' + sP.plotPath + 'job_fof.bsub'
    endif
  
  endfor ;snapRange

end

; makeQueueTestJobs():

pro makeQueueTestJobs

  ; config
  spawnJobs = 1
  nProcs    = 8 ; per node
  queueName = 'nancy'
  jobPath   = '/n/home07/dnelson/queue.test/'
  
  ; get host group name
  spawn, 'bqueues -l '+queueName, result
  pos_start = strpos(result,'HOSTS:')
  
  w = ( where(pos_start eq 0,count) )[0]
  if count eq 0 then message,'Error'
  pos_end = strpos(result[w],'/')
  
  hostGroupName = strmid(result[w],strlen('HOSTS:'),strlen(result[w])-strlen('HOSTS:')-2)
  hostGroupName = str( hostGroupName )
  
  ; get list of nodes in this host group
  maxNodeNameLen = 10
  maxHosts = 500
  
  spawn,'bhosts ' + hostGroupName, result
  result = result[1:*]
  
  hosts = strarr(n_elements(result) < maxHosts)
  
  for i=0,n_elements(hosts)-1 do $
    hosts[i] = str( strmid(result[i],0,maxNodeNameLen) )
  
  print,'Using host group ['+hostGroupName+'], running on ['+str(n_elements(hosts))+'] hosts:'
  print,hosts
  print,''
  
  ; loop over all hosts, submit test job  
  foreach hostName,hosts,i do begin
  
    ; write bjob file
    jobFileName = jobPath + 'job_test.bsub'
    
    print,'['+str(i)+'] Writing: '+hostName
    
    ; host exclusion string
    selectStr = ''
    if n_elements(excludeHosts) gt 0 then begin
      selectStr = 'select['
      for i=0,n_elements(excludeHosts)-1 do selectStr += 'hname!='+excludeHosts[i]+' && '
      selectStr = strmid(selectStr,0,strlen(selectStr)-4) + '] '
    endif
      
    openW, lun, jobFilename, /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#BSUB -q nancy'
    printf,lun,'#BSUB -J ' + hostName
    printf,lun,'#BSUB -n ' + str(nProcs)
    printf,lun,'#BSUB -m ' + str(hostName)
    printf,lun,'#BSUB -x'
    printf,lun,'#BSUB -g /dnelson/orbit'
    printf,lun,'#BSUB -o run.txt'
    printf,lun,'#BSUB -e run.err'
    printf,lun,'#BSUB -cwd ' + jobPath
    printf,lun,''
    printf,lun,'python test_script.py >> ${LSB_JOBNAME}.txt'
    
    ; close
    close,lun
    free_lun,lun
      
    ; add to queue if requested
    if (spawnJobs) then begin
      spawn, 'bsub < ' + jobPath + 'job_test.bsub', result
      print,'  '+result
      
      ; file cleanup
      wait,0.2
      spawn, 'rm ' + jobPath + 'job_test.bsub'
    endif
    
  endforeach

end

; plotQueueTest()

pro plotQueueTest

  ; config
  jobPath   = '/n/home07/dnelson/queue.test/'
  
  host     = []
  testtime = []
  datetime = []
  
  files = file_search(jobPath + 'hero*.txt')
  
  ; load
  foreach file,files do begin
    nRows = file_lines(file)
    if nRows eq 0 then continue
    
    ; hostname
    hostname = strsplit( file, '/', /extract)
    hostname = hostname[-1]
    hostname = strsplit( hostname, '.txt', /extract)
    hostname = hostname[0]
    
    openR, lun, file, /GET_LUN
    
    ; add each row to keeper arrays
    for i=0,nRows-1 do begin
      row = ''
      readF, lun, row
        
      ; parse
      row = strsplit(row, ", '", /extract)
      test = float( strmid(row[0],1) )
      time = row[1] + ' ' + row[2]
      
      host     = [host, hostname]
      testtime = [testtime, test]
      datetime = [datetime, time]
    endfor
      
    close, lun
    free_lun, lun
    
  endforeach
  
  print,'Found ['+str(n_elements(host))+'] entries.'
  
  start_PS, jobPath + 'histo.eps'
    plothist, testTime / min(testTime), /auto, xtitle=textoidl("\Delta t / min ( \Delta t )"), $
      ytitle="Number of Nodes", /ylog
  end_PS
  
  ; sort
  sort_inds = reverse( sort(testTime) )
  testTime = testTime[sort_inds]
  dateTime = dateTime[sort_inds]
  host     = host[sort_inds]
  
  ; print 10 worst offenders
  print,'Top ten worst times:'
  for i=0,9 do print,str(testTime[i]) + '  ('+host[i]+')'

  stop
  
end