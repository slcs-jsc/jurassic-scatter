#include "jurassic.h"
#include "control.h"
#include "misc.h"

int main(int argc, char *argv[]) {
  
  static ctl_t ctl;
  static tbl_t *tbl;
  
  FILE *out;
  
  char filename[2*LEN];
  
  int id, ig, ip, it;
  
  /* Check arguments... */
  if(argc<2)
    ERRMSG("Give parameters: <ctl>");
  
  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Allocate... */
  ALLOC(tbl, tbl_t, 1);
  
  /* Read tables... */
  init_tbl(&ctl, tbl);
  
  /* Loop over emitters... */
  for(ig=0; ig<ctl.ng; ig++) {

    /* Loop over channels... */
    for(id=0; id<ctl.nd; id++) {
      
      /* Set filename... */
      sprintf(filename, "%s_%.4f_%s.bin",
	      ctl.tblbase, ctl.nu[id], ctl.emitter[ig]);
      
      /* Create file... */
      LOGMSG(2, printf("Write binary table: %s\n", filename));
      if(!(out=fopen(filename, "w")))
        ERRMSG("Cannot create file!");
      
      /* Write data... */
      FWRITE(&tbl->np[ig][id], int, 1, out);
      FWRITE(&tbl->p[ig][id], double, tbl->np[ig][id], out);
      FWRITE(tbl->nt[ig][id], int, tbl->np[ig][id], out);
      for(ip=0; ip<tbl->np[ig][id]; ip++) {
	FWRITE(&tbl->t[ig][id][ip], double, tbl->nt[ig][id][ip], out);
	FWRITE(tbl->nu[ig][id][ip], int, tbl->nt[ig][id][ip], out);
	for(it=0; it<tbl->nt[ig][id][ip]; it++) {
	  FWRITE(&tbl->u[ig][id][ip][it], float,
		 tbl->nu[ig][id][ip][it], out);
	  FWRITE(&tbl->eps[ig][id][ip][it], float,
		 tbl->nu[ig][id][ip][it], out);
	}
      }
      
      /* Close file... */
      fclose(out);
    }
  }
  
  /* Free... */
  free(tbl);
  
  return EXIT_SUCCESS;
}
