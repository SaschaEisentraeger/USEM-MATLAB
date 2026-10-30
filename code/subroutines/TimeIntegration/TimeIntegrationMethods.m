function methods = TimeIntegrationMethods()
%TIMEINTEGRATIONMETHODS Table of the direct time integration methods.
%
%   methods = TimeIntegrationMethods() returns a struct array with the
%   fields
%     code    method code (settings.TimeIntegration.Method)
%     name    description
%     handle  function that integrates the equations of motion
%     args    function handle opt -> cell array of the method-specific
%             trailing arguments; opt is settings.TimeIntegration with the
%             additional field ElemData (element matrices, method 109)
%
%   All methods share the calling convention
%
%     [T,U] = handle(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,ActiveDoFs,...
%                    ExciteParam,SignalType,dtmax,tsim,args(opt){:})
%
%   K, M, C       system matrices of the active DOFs
%   F, F_dt       constant load and time-dependent load parts (struct with
%                 fields Neumann, K, C, M, multiplied by s, s, ds/dt and
%                 d2s/dt2 of the excitation signal s)
%   U0, V0        initial displacements and velocities
%   NumDoF        number of all DOFs
%   DOFsave       DOFs (of the reduced system) whose history is returned
%   ActiveDoFs    active DOFs; if all of them are saved, U is expanded to
%                 all NumDoF DOFs
%   dtmax         time step candidates [dtf, dte, dtImp]: excitation,
%                 stability limit of the central difference method, and
%                 prescribed step of the implicit methods
%   tsim          simulation time
%   T, U          time instants and displacement history
%
%   To add a method, write a function with this signature (see
%   NewmarkMethod for a template that uses the shared helpers in
%   TimeIntegration/Common) and add one row to the table below.
%   Solver_Transient looks the method up by its code.

none = @(opt) {};

table = {
% code  name                                          handle                            trailing arguments
%-------------------------------------------------------------------------- explicit, 2nd order
  101, 'Central difference method',                   @CentralDifferenceMethod,         @(opt) {true}
  102, 'Noh-Bathe method',                            @NohBatheMethod,                  none
  103, 'Soares method (explicit)',                    @SoaresMethodExplicit,            none
  104, 'Mirbagheri method',                           @MirbagheriMethod,                none
  105, 'Kim method (explicit)',                       @KimMethodExplicit,               none
  106, 'Zhang-Xing method 1',                         @ZhangXingMethod1,                none
  107, 'Zhang-Xing method 2',                         @ZhangXingMethod2,                none
  108, 'Explicit generalized-alpha method',           @ExplicitGeneralizedAlphaMethod,  none
  109, 'Adaptive hybrid explicit method (Soares)',    @AdaptiveHybridExplicitMethod,    @(opt) {opt.ElemData}
%-------------------------------------------------------------------------- explicit, higher order
  201, 'Kim method (explicit, order 3 or 4)',         @KimExplicit_ho,                  @(opt) {opt.Order}
  202, 'Two-stage Kim method (explicit)',             @KimMethod2Stage,                 none
%-------------------------------------------------------------------------- implicit, 2nd order
  301, 'Newmark method',                              @NewmarkMethod,                   none
  302, 'HHT-alpha method',                            @HHTAlphaMethod,                  none
  303, 'Generalized-alpha method',                    @GeneralizedAlphaMethod,          none
  304, 'Bathe method (Bathe-Baig)',                   @BatheBaigMethod,                 none
  305, 'beta1-beta2-Bathe method',                    @b1b2_BatheMethod,                none
  306, 'rho-Bathe method',                            @rho_BatheMethod,                 none
  307, 'Soares method (implicit)',                    @SoaresMethodImplicit,            none
  308, 'Kim method (implicit)',                       @KimMethodImplicit,               none
  309, 'Kim method (implicit, collocation)',          @KimMethodImplicitCollocation,    none
%-------------------------------------------------------------------------- implicit, higher order
  401, 'Kim method (Lagrange, order 1 to 5)',         @KimLagrange_ho,                  @(opt) {opt.Order}
  402, 'Soares method (implicit, higher order)',      @SoaresMethodImplicit_ho,         none
  403, 'Pade-based method (order 2 to 10)',           @PadeMethod,                      @(opt) {opt.Order,opt.rho,opt.nInterp}
%-------------------------------------------------------------------------- other
  501, 'Highly precise time integration (HPTIM)',     @HPTIM_Newmark,                   @(opt) {opt.Order}
  601, 'Central difference method (parallel code)',   @CentralDifferenceMethod,         @(opt) {true}
};

methods = cell2struct(table,{'code','name','handle','args'},2);

end
