function trajinfo(psffile, dcdfile)

      ta = mdload(psffile);
      ta, boxsize = mdload(dcdfile, top=ta);

      return ta, boxsize 

end 
