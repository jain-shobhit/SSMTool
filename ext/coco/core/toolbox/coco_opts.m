function coco_opts()
%List of options used by the continuation core.
%
%The options structure has a depth of three. The first level consists of
%'classes'. The second level consists of 'properties' belonging to the
%corresponding 'classes'. The third level consists of 'values' of the
%corresponding 'properties'.
%
%The following classes and properties are recognised by COCO (defaults in
%square brackets):
%
%Class pdat:
%   The class pdat (problem data) stores information about the current
%   problem provided by the user.
%
%   F            : function F(x,p) as provided by the user
%   x0           : initial value for x
%   par          : initial value for p
%   cont_par     : index of continuation parameter
%   cont_par_int : interval for continuation parameter
%   vectorised   : user function F is vectorised [ 'true' ]
% 
%Class func:
%   The class func represents the equation for which families of solutions
%   will be computed. This is typically an algorithm in functional form.
%
%   F        : algorithm in functional form, this property is used,
%              but not set by COCO
%   DFDX     : linearisation of algorithm, this property is used,
%              but not set by COCO
%   DFDP     : derivative of algorithm with respect to parameters, this
%              property is used, but not set by COCO
%   linsolve : solver that solves linear equations of the form
%              DFDX(x,p)*d = F(x,p) [ @coco_default_linsolve ]
% 
%Class xfunc:
%   The class xfunc represents an extended system that embedds the
%   algorithm defined in class func.
%
%   F        : extended system embedding the algorithm in func [ @xfunc_F ]
%   DFDX     : linearisation of extended system [ @xfunc_DFDX ]
%   linsolve : solver that solves extended linear equations of the form
%              DFDX(x,p)*d = F(x,p) [ @coco_default_linsolve ]
%   update   : update last-point information [ @xfunc_update ]
%              Note: this property is not yet implemented.
% 
%Class nwtn:
%   The class nwtn represents the Newton corrector. Its properties store
%   variables that influence the behaviour of Newton's method, as well as
%   intermediate computational results.
%
%* Properties that influence the algorithm's behaviour:
%
%   ItMX     [     10 ] : max. number of iterations
%   ItNW     [     [] ] : max. number of full Newton iterations
%   SubItMX  [      8 ] : max. number of damping steps
%   TOL      [ 1.0e-8 ] : convergence criterion
%   LogLevel [      1 ] : level of diagnostic output
%   MaxStep  [    0.1 ] : max. size of Newton step
%   ga0      [    1.0 ] : initial damping factor
%   al       [    0.5 ] : decrease damping factor
%
%* Properties defining functions callback handlers:
%
%   print_headline      [ @coco_default_print    ]
%   print_data          [ @coco_default_print    ]
%   begin_callback      [ @coco_default_callback ]
%   begin_step_callback [ @coco_default_callback ]
%   end_callback        [ @coco_default_callback ]
%   end_step_callback   [ @coco_default_callback ]
% 
%* Property defining the equations to be solved:
%
%   func [ 'func' ] : functions to use
%
%* Properties used for internal data:
%
%   It         : number of completed Newton iterations
%   accept     : flag if iteration converged
%   ftm        : total time spent on evaluating F
%   dftm       : total time spent on evaluating DFDX
%   stm        : total time spent solving the linearised equation
%   SubIt      : number of completed damping steps
%   x0         : initial guess for solution of F(x,p)=0
%   x          : current solution
%   f          : current residual
%   norm_f_old : norm of previously computed residual
%   J          : current value of DFDX
%   d          : current correction vector
%   ga         : current damping factor
%
%Class cont:
%   The class cont represents the generic continuation method. Its
%   properties store variables that influence its behaviour, as well as
%   intermediate computational results.
%
%
%* Properties that influence the algorithm's behaviour:
%
%   MaxRes    [   0.1 ] : max. residuum for prediction step
%   al_max    [   7.0 ] : max. angle between two consecutive tangents
%   h0        [   0.1 ] : initial continuation step size
%   h_max     [   0.5 ] : max. continuation step size
%   h_min     [  0.01 ] : min. continuation step size
%   h_fac_min [   0.5 ] : min. step size adaption factor
%   h_fac_max [   2.0 ] : max. step size adaption factor
%   ga        [  0.95 ] : adaption security factor
%   ItMX      [   100 ] : max. number of continuation steps
%   LogLevel  [ [1 0] ] : diagnostic output level
%   NPR       [    10 ] : primnt point information every NPR steps
%   NSV       [    [] ] : save solution every NSV steps, use NPR if empty
%
%* Properties defining functions and callback handlers:
%
%   continuer      [              'state' ] class name of finite state
%                                           machine for continuation
%   corrector      [ @coco_nwtn_step      ]
%   print_headline [ @coco_default_print  ]
%   print_data     [ @coco_default_print  ]
%   save_full      [ @coco_default_save   ]
%   save_reduced   [ @coco_default_save   ]
%   update         [ @coco_default_update ]
% 
%* Property defining the equations to be solved:
%
%   func  [  'func' ] algorithm for which families solutions are computed
%   xfunc [ 'xfunc' ] extended function embedding the algorithm in func
%
%* Properties used for internal data:
%
%   state     : current state of finite-state machine
%   accept    : flag if continuation terminated with success
%   branch    : flag for direction from initial point
%   xidx      : indices of solution components in the extended solution
%               vector
%   pidx      : indices of continuation parameters in extended solution
%               vector
%   It        : number of continuation steps
%   h         : current continuation step size
%   lab       : solution label
%   next_lab  : next solution label that can be allocated
%   tm        : start time of continuation
%   x0        : initial solution
%   p0        : initial parameters
%   cont_par  : index of continuation parameter
%   pint      : parameter interval
%   u         : extended solution vector
%   u0        : initial extended solution vector
%   p         : current parameter vector
%   us        : vector tangent to family in u
%   us0       : vector tangent to family in u0
%   v         : predicted extended solution vector
%   spt       : special point type
% 
%Class sol:
%   The class sol is used to store the current solution and associated
%   data.
%
%   p   : current parameter values
%   u   : current extended solution vector
%   us  : tangent vector at u
%   spt : special point type
%   
%Class state:
%   The class state represents the finite state machine that encodes the
%   continuation method to be used. Its properties are states, which are
%   responsible for performing specific tasks of a continuation algorithm.
%
%   init            [ @state_init            ]
%   predict         [ @state_predict         ]
%   correct         [ @state_correct         ]
%   check_solution  [ @state_check_solution  ]
%   compute_events  [ @state_compute_events  ]
%   handle_events   [ @state_handle_events   ]
%   insert_point    [ @state_insert_point    ]
%   update          [ @state_update          ]
%   refine_step     [ @state_refine_step     ]
%   handle_boundary [ @state_handle_boundary ]
% 
%Class pred:
%   The class pred represents a specific prediction method.
%
%   init    [ @tpred_init    ]
%   update  [ @tpred_update  ]
%   predict [ @tpred_predict ]
% 
%Class bddat:
%   The class bddat provides functions for manipulating the data stored in
%   a bifurcation diagram.
%
%   init     [ @bddat_init    ]
%   insert   [ @bddat_insert  ]
%   append   [ @bddat_append  ]
%   prepend  [ @bddat_prepend ]
%   usr_init                    : not set by COCO
%   usr_data                    : not set by COCO
% 
%Class run:
%   The class run collects some useful information about the current run.
%
%   run       : label for run
%   mfname    : filename containing user function
%   dir       : directory for saving data
%   bdfname   : name of file for saving bifurcation diagram
%   isol_type : type of initial solution
%   sol_type  : type of sought solution
% 
%Class toolbox:
%   The class toolbox collects information about a toolbox. The data of
%   this class is saved in all solution files for identification purposes.
%   With this information a specific toolbox should be able to decide if it
%   can handle the data in a solution file. Toolboxes should add small
%   items of data if needed for identification.
%
%   name : name of the toolbox called by COCO. This is set by COCO.
% 
%Class ptlist:
%   The class ptlist is used for temporarily storing a list of special
%   points that were found in a continuation step. The points are stored so
%   they can be ordered before processing. Ptlist is an array of structs
%   with the following elements:
%
%   p   : current parameter values
%   u   : current extended solution vector
%   us  : tangent vector at u
%   spt : special point type
% 
%NOTE: Any continuation toolbox compatible with COCO should provide
%information about additional options used by this toolbox. Type 'help
%TOOLBOX_opts', where TOOLBOX is the name of a continuation toolbox.
%

help coco_opts

