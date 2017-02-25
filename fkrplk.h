      parameter(nemax=30,nmumax=15,ntaumax=30)
      common /univ/ phi(0:nemax,-nmumax:nmumax,0:ntaumax,0:1)
      common /univ1/ amu(-nmumax:nmumax), ne, nmu, ntau, 
     +       beta(0:nemax), rlmbda(0:ntaumax), dy