// STATUS BAR  
// Modified
// http://www.rosshemsley.co.uk/2011/02/creating-a-progress-bar-in-c-or-any-other-console-app/

static inline void loadBar(int x, int n, int r, int w)
{
  // Only update r times.
  if ( x % (n/r) != 0 ) return;

  // Calculuate the correlation of complete-to-incomplete.
  double correlation = x/(double)n;
  int c = correlation * w;

  // Show the percentage complete.
  printf("%3d%%[", (int)(correlation*100) );

  // Show the load bar.
  for (int x=0; x<c; x++)
    printf("=");

  for (int x=c; x<w; x++)
    printf(" ");

  // ANSI Control codes to go back to the
  // previous line and clear it.
  // printf("]\n33[F33[J");

  printf("\r"); // Move to the first column
  fflush(stdout);
}
