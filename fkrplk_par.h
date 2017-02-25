      parameter(nemax=60,nmumax=30,ntaumax=30)
      common /univ/ phi(0:1,0:nemax,-nmumax:nmumax,0:ntaumax)
      common /univ1/ amu(-nmumax:nmumax), ne, nmu, ntau, 
     +       beta(0:nemax), rlmbda(0:ntaumax), dy