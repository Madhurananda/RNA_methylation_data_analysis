import os,sys
import glob

if __name__=='__main__':
   if len(sys.argv)<3:
      print("Usage: {} {} {}".format(sys.argv[0], "basefolder", "basefile"))
      print("Example:")
      print("       {} {} {}".format(sys.argv[0], "Curlcake_m6a", "1204665"))
      print("       {} {} {}".format(sys.argv[0], "Curlcake_m6a", "1368282"))
      print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1204670"))
      print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1368288"))
      sys.exit(1)

   basefolder = sys.argv[1]
   basetar = sys.argv[2]
   
   print('basefolder:', basefolder)
   print('basetar: ', basetar)

   if len(sys.argv)>3:
#       print('\n**** Just checking ****\n')
      start_fid = int(sys.argv[3])
      end_fid = int(sys.argv[4])
      files = []
      for _i in range(start_fid, end_fid):
         if os.path.isfile( basefolder+"/"+basetar+"-"+str(_i)+".fast5.tar"):
            files.append( basefolder+"/"+basetar+"-"+str(_i)+".fast5.tar")
   else:
      print('There are exactly 3 inputs.')
      print('command: ', os.path.join(basefolder, "*.tar.*"))
      files = glob.glob(os.path.join(basefolder, "*.tar.*"))
#       print('command: ', os.path.join(basefolder, basetar+"*.tar.*"))
#       files = glob.glob(os.path.join(basefolder, basetar+"*.tar.*"))
#       files = glob.glob(os.path.join(basefolder, basetar))
#       files = glob.glob(os.path.join(basefolder, basetar+"*.tar"))

   print('files: ', files)

   if not os.path.isdir(basefolder+'/'+basetar):
      os.system("mkdir {}".format( basefolder+'/'+basetar));

   for _f in files:
      _fid = _f.split('-')[-1].split('.')[0]
      print( '_fid: ', _fid )
#       if not os.path.isdir( "{}".format( basefolder+'/'+basetar+'/'+_fid ):
#          os.system("mkdir {}".format( basefolder+'/'+basetar+'/'+_fid));
#       print("tar -xf {} -C {}".format( _f, basefolder+'/'+basetar+'/'+_fid ))
#       os.system( "tar -xf {} -C {}".format( _f, basefolder+'/'+basetar+'/'+_fid ) )
#       os.system( "/usr/bin/rm {}".format( _f))

