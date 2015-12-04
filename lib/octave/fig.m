function fig (name,ver,plotting)

  ## usage:  fig (name,plotting)
  ##
  ## Generates name_ver.eps if plotting==1
  ## Does nothing if plotting==0

if nargin>1
  name = sprintf("%s_%s",name,ver);
endif

if nargin<3
  plotting = 1;
endif

# if nargin>1
#   if isreal(ver)
#     ver = "";
#     plotting = ver;
#   endif
# endif


if (plotting>0)
  if plotting==1
    colour="-mono";
  else
    colour="-color";
  endif
  cd Figs
  ##figname = sprintf("%s.fig",name)
 figname = sprintf("%s.eps",name);
## set(gca(),"fontsize",30) 
 print(figname,colour, "-deps");
  ##print(figname,colour, "-F:14");
  ##system(sprintf("(fig2eps --bbox=dvips --forcespecial --nogv %s.fig) &>/dev/null",name));
  ##system(sprintf("eps2eps  %s.eps junk.eps",name));
  ##system(sprintf("mv junk.eps %s.eps",name));
  ##system(sprintf("rm -f %s.fig",name));
  cd ..
else
  printf("Not saving figure %s\n",name);
endif

endfunction