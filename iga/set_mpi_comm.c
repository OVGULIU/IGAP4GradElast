void set_mpi_comm(
int *nselfIdx_,
int **selfIdx_,
int *nsendPart_,
int **sendPart_,
int *nsendIdx_,
int **sendIdx_,
int **sendPtr_,
int *nrecvPart_,
int **recvPart_,
int *nrecvIdx_,
int **recvIdx_,
int **recvPtr_,
int ndof,
int porder,
int rank,
int npart[],
int ipart[],
int nelem[],
int nbasis[]
)
{
    // ---- 
    int nselfIdx;
    int nsendPart;
    int nsendIdx;
    int nrecvPart;
    int nrecvIdx;

    // ---- self
    nselfIdx=ndof*(nbasis[0]*nbasis[1]*nbasis[2]);

    // ---- send
    nsendPart=(int)(ipart[2] != 0)
             +(int)(ipart[1] != 0)
             +(int)(ipart[0] != 0)
             +(int)(ipart[1] != 0 && ipart[2] != 0)
             +(int)(ipart[0] != 0 && ipart[2] != 0)
             +(int)(ipart[0] != 0 && ipart[1] != 0)
             +(int)(ipart[0] != 0 && ipart[1] != 0 && ipart[2] != 0);
    nsendIdx=0;
    if (ipart[2] != 0)                                   { nsendIdx+=ndof*(nbasis[0]*nbasis[1]*porder); }
    if (ipart[1] != 0)                                   { nsendIdx+=ndof*(nbasis[0]*nbasis[2]*porder); }
    if (ipart[0] != 0)                                   { nsendIdx+=ndof*(nbasis[1]*nbasis[2]*porder); }
    if (ipart[1] != 0 && ipart[2] != 0)                  { nsendIdx+=ndof*(nbasis[0]*porder*porder); }
    if (ipart[0] != 0 && ipart[2] != 0)                  { nsendIdx+=ndof*(nbasis[1]*porder*porder); }
    if (ipart[0] != 0 && ipart[1] != 0)                  { nsendIdx+=ndof*(nbasis[2]*porder*porder); }
    if (ipart[0] != 0 && ipart[1] != 0 && ipart[2] != 0) { nsendIdx+=ndof*(porder*porder*porder); }

    // ---- recv
    nrecvPart=(int)(ipart[2] != npart[2]-1)
             +(int)(ipart[1] != npart[1]-1)
             +(int)(ipart[0] != npart[0]-1)
             +(int)(ipart[1] != npart[1]-1 && ipart[2] != npart[2]-1)
             +(int)(ipart[0] != npart[0]-1 && ipart[2] != npart[2]-1)
             +(int)(ipart[0] != npart[0]-1 && ipart[1] != npart[1]-1)
             +(int)(ipart[0] != npart[0]-1 && ipart[1] != npart[1]-1 && ipart[2] != npart[2]-1);
    nrecvIdx=0;
    if (ipart[2] != npart[2]-1)                                                     { nrecvIdx+=ndof*(nbasis[0]*nbasis[1]*porder); }
    if (ipart[1] != npart[1]-1)                                                     { nrecvIdx+=ndof*(nbasis[0]*nbasis[2]*porder); }
    if (ipart[0] != npart[0]-1)                                                     { nrecvIdx+=ndof*(nbasis[1]*nbasis[2]*porder); }
    if (ipart[1] != npart[1]-1 && ipart[2] != npart[2]-1)                           { nrecvIdx+=ndof*(nbasis[0]*porder*porder); }
    if (ipart[0] != npart[0]-1 && ipart[2] != npart[2]-1)                           { nrecvIdx+=ndof*(nbasis[1]*porder*porder); }
    if (ipart[0] != npart[0]-1 && ipart[1] != npart[1]-1)                           { nrecvIdx+=ndof*(nbasis[2]*porder*porder); }
    if (ipart[0] != npart[0]-1 && ipart[1] != npart[1]-1 && ipart[2] != npart[2]-1) { nrecvIdx+=ndof*(porder*porder*porder); }


    // ---- 
    int *selfIdx  = (int*) malloc(nselfIdx*sizeof(int));
    int *sendPart = (int*) malloc(nsendPart*sizeof(int));
    int *sendIdx  = (int*) malloc(nsendIdx*sizeof(int));
    int *sendPtr  = (int*) malloc((nsendPart+1)*sizeof(int));
    int *recvPart = (int*) malloc(nrecvPart*sizeof(int));
    int *recvIdx  = (int*) malloc(nrecvIdx*sizeof(int));
    int *recvPtr  = (int*) malloc((nrecvPart+1)*sizeof(int));

    int ibasis,ia;
    //int nbasis_x_active=nelem[0]+porder;
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    
    // ---- selfIdx
    ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { selfIdx[ia++]=ndof*ibasis+idof; }
    }}}

    // ---- sendPart & sendPtr
    nsendIdx=0;
    ia=0;
    if (ipart[2] != 0)                                   {sendPart[ia]=rank-1                           ; sendPtr[ia]=nsendIdx; nsendIdx+=ndof*(nbasis[0]*nbasis[1]*porder);ia++;}
    if (ipart[1] != 0)                                   {sendPart[ia]=rank-npart[2]                    ; sendPtr[ia]=nsendIdx; nsendIdx+=ndof*(nbasis[0]*nbasis[2]*porder);ia++;}
    if (ipart[0] != 0)                                   {sendPart[ia]=rank-npart[1]*npart[2]           ; sendPtr[ia]=nsendIdx; nsendIdx+=ndof*(nbasis[1]*nbasis[2]*porder);ia++;}
    if (ipart[1] != 0 && ipart[2] != 0)                  {sendPart[ia]=rank-npart[2]-1                  ; sendPtr[ia]=nsendIdx; nsendIdx+=ndof*(nbasis[0]*porder*porder);ia++;}
    if (ipart[0] != 0 && ipart[2] != 0)                  {sendPart[ia]=rank-npart[1]*npart[2]-1         ; sendPtr[ia]=nsendIdx; nsendIdx+=ndof*(nbasis[1]*porder*porder);ia++;}
    if (ipart[0] != 0 && ipart[1] != 0)                  {sendPart[ia]=rank-npart[1]*npart[2]-npart[2]  ; sendPtr[ia]=nsendIdx; nsendIdx+=ndof*(nbasis[2]*porder*porder);ia++;}
    if (ipart[0] != 0 && ipart[1] != 0 && ipart[2] != 0) {sendPart[ia]=rank-npart[1]*npart[2]-npart[2]-1; sendPtr[ia]=nsendIdx; nsendIdx+=ndof*(porder*porder*porder);ia++;}
    sendPtr[ia]=nsendIdx;

    // ---- sendIdx
    ia=0;
    // 0 (face0)
    if (ipart[2] != 0){
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 1 (face1)
    if (ipart[1] != 0){
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 2 (face2)
    if (ipart[0] != 0){
    for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 3 (edge0)
    if (ipart[1] != 0 && ipart[2] != 0){
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 4 (edge1)
    if (ipart[0] != 0 && ipart[2] != 0){
    for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 5 (edge2)
    if (ipart[0] != 0 && ipart[1] != 0){
    for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 6 (corner)
    if (ipart[0] != 0 && ipart[1] != 0 && ipart[2] != 0){
    for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
    }}}}

    // ---- recvPart & recvPtr
    nrecvIdx=0;
    ia=0;
    if (ipart[2] != npart[2]-1)                                                     {recvPart[ia]=rank+1                           ; recvPtr[ia]=nrecvIdx; nrecvIdx+=ndof*(nbasis[0]*nbasis[1]*porder);ia++;}
    if (ipart[1] != npart[1]-1)                                                     {recvPart[ia]=rank+npart[2]                    ; recvPtr[ia]=nrecvIdx; nrecvIdx+=ndof*(nbasis[0]*nbasis[2]*porder);ia++;}
    if (ipart[0] != npart[0]-1)                                                     {recvPart[ia]=rank+npart[1]*npart[2]           ; recvPtr[ia]=nrecvIdx; nrecvIdx+=ndof*(nbasis[1]*nbasis[2]*porder);ia++;}
    if (ipart[1] != npart[1]-1 && ipart[2] != npart[2]-1)                           {recvPart[ia]=rank+npart[2]+1                  ; recvPtr[ia]=nrecvIdx; nrecvIdx+=ndof*(nbasis[0]*porder*porder);ia++;}
    if (ipart[0] != npart[0]-1 && ipart[2] != npart[2]-1)                           {recvPart[ia]=rank+npart[1]*npart[2]+1         ; recvPtr[ia]=nrecvIdx; nrecvIdx+=ndof*(nbasis[1]*porder*porder);ia++;}
    if (ipart[0] != npart[0]-1 && ipart[1] != npart[1]-1)                           {recvPart[ia]=rank+npart[1]*npart[2]+npart[2]  ; recvPtr[ia]=nrecvIdx; nrecvIdx+=ndof*(nbasis[2]*porder*porder);ia++;}
    if (ipart[0] != npart[0]-1 && ipart[1] != npart[1]-1 && ipart[2] != npart[2]-1) {recvPart[ia]=rank+npart[1]*npart[2]+npart[2]+1; recvPtr[ia]=nrecvIdx; nrecvIdx+=ndof*(porder*porder*porder);ia++;}
    recvPtr[ia]=nrecvIdx;
    
    // ---- recvIdx
    ia=0;
    // 0 (face0)
    if (ipart[2] != npart[2]-1){
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 1 (face1)
    if (ipart[1] != npart[1]-1){
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 2 (face2)
    if (ipart[0] != npart[0]-1){
    for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 3 (edge0)
    if (ipart[1] != npart[1]-1 && ipart[2] != npart[2]-1){
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
    for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 4 (edge1)
    if (ipart[0] != npart[0]-1 && ipart[2] != npart[2]-1){
    for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 5 (edge2)
    if (ipart[0] != npart[0]-1 && ipart[1] != npart[1]-1){
    for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
    for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
    }}}}
    // 6 (corner)
    if (ipart[0] != npart[0]-1 && ipart[1] != npart[1]-1 && ipart[2] != npart[2]-1){
    for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
    for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
    for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
    }}}}


    // ----
    *nselfIdx_      = nselfIdx;
    *nsendPart_     = nsendPart;
    *nsendIdx_      = nsendIdx;
    *nrecvPart_     = nrecvPart;
    *nrecvIdx_      = nrecvIdx;

    *selfIdx_       = selfIdx;
    *sendPart_      = sendPart;
    *sendIdx_       = sendIdx;
    *sendPtr_       = sendPtr;
    *recvPart_      = recvPart;
    *recvIdx_       = recvIdx;
    *recvPtr_       = recvPtr;
}
