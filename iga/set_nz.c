void set_nz(int *d_nz, int *o_nz, int nbasis[], int npart[], int ipart[], int ndof, int porder)
{
    int ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
    for (int idof=0; idof<ndof; idof++) {
        d_nz[ia]=ndof*( ( ibasis_x < porder ? ibasis_x : porder )+1+( ibasis_x > (nbasis[0]-1)-porder ? (nbasis[0]-1)-ibasis_x:porder) )
                     *( ( ibasis_y < porder ? ibasis_y : porder )+1+( ibasis_y > (nbasis[1]-1)-porder ? (nbasis[1]-1)-ibasis_y:porder) )
                     *( ( ibasis_z < porder ? ibasis_z : porder )+1+( ibasis_z > (nbasis[2]-1)-porder ? (nbasis[2]-1)-ibasis_z:porder) );
        o_nz[ia]=ndof*( ( (ipart[0]==0 && ibasis_x < porder) ? ibasis_x : porder )+1+( (ipart[0]==npart[0]-1 && ibasis_x > (nbasis[0]-1)-porder ) ? (nbasis[0]-1)-ibasis_x : porder ) )
                     *( ( (ipart[1]==0 && ibasis_y < porder) ? ibasis_y : porder )+1+( (ipart[1]==npart[1]-1 && ibasis_y > (nbasis[1]-1)-porder ) ? (nbasis[1]-1)-ibasis_y : porder ) )
                     *( ( (ipart[2]==0 && ibasis_z < porder) ? ibasis_z : porder )+1+( (ipart[2]==npart[2]-1 && ibasis_z > (nbasis[2]-1)-porder ) ? (nbasis[2]-1)-ibasis_z : porder ) )
                 -d_nz[ia];
        ia++;
    }}}}
}
