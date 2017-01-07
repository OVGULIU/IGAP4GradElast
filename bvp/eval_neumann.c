void eval_neumann(double *traction, double par_neumann[], int face_id, int bc_order, int idof, double xq[])
{
    int ndof=3;

    //                                          0                          |                          1                               // bc_order
    //                -------------------------------------------------------------------------------------------------------------------------
    //                        0        ,        1        ,        2        |        0        ,        1        ,        2             // idof
    //                -------------------------------------------------------------------------------------------------------------------------
    double neumann[]={        0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 0
                              0        ,  par_neumann[0] ,  par_neumann[0] ,        0        ,        0        ,        0        ,    // face 1
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 2
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 3
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 4
                              0        ,        0        ,        0        ,        0        ,        0        ,        0         };  // face 5
    
    *traction=neumann[(2*ndof)*face_id+ndof*bc_order+idof];
}
