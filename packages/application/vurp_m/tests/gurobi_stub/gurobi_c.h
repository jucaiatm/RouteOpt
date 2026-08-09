#pragma once

// Minimal declarations used only by the VURP-M syntax-only CI target.
// Production builds must use the official Gurobi header and library.

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GRBenv GRBenv;
typedef struct GRBmodel GRBmodel;

#define GRB_INFINITY 1.0e100

#define GRB_LESS_EQUAL '<'
#define GRB_GREATER_EQUAL '>'
#define GRB_EQUAL '='
#define GRB_CONTINUOUS 'C'

#define GRB_MINIMIZE 1
#define GRB_METHOD_DUAL 1

#define GRB_OPTIMAL 2
#define GRB_INFEASIBLE 3
#define GRB_INF_OR_UNBD 4

#define GRB_INT_ATTR_NUMVARS "NumVars"
#define GRB_INT_ATTR_NUMCONSTRS "NumConstrs"
#define GRB_INT_ATTR_MODELSENSE "ModelSense"
#define GRB_INT_ATTR_STATUS "Status"

#define GRB_DBL_ATTR_OBJVAL "ObjVal"
#define GRB_DBL_ATTR_X "X"
#define GRB_DBL_ATTR_PI "Pi"
#define GRB_DBL_ATTR_LB "LB"
#define GRB_DBL_ATTR_UB "UB"
#define GRB_DBL_ATTR_OBJ "Obj"

#define GRB_INT_PAR_OUTPUTFLAG "OutputFlag"
#define GRB_INT_PAR_METHOD "Method"
#define GRB_INT_PAR_THREADS "Threads"
#define GRB_INT_PAR_NONCONVEX "NonConvex"

const char *GRBgeterrormsg(GRBenv *env);
GRBenv *GRBgetenv(GRBmodel *model);

int GRBloadenv(GRBenv **env, const char *logfilename);
void GRBfreeenv(GRBenv *env);

int GRBnewmodel(GRBenv *env, GRBmodel **model, const char *name, int numvars,
                double *obj, double *lb, double *ub, char *vtype,
                char **varnames);
void GRBfreemodel(GRBmodel *model);

int GRBupdatemodel(GRBmodel *model);
int GRBoptimize(GRBmodel *model);

int GRBsetintparam(GRBenv *env, const char *paramname, int value);
int GRBsetintattr(GRBmodel *model, const char *attrname, int value);
int GRBsetdblattrelement(GRBmodel *model, const char *attrname,
                         int element, double value);

int GRBgetintattr(GRBmodel *model, const char *attrname, int *value);
int GRBgetdblattr(GRBmodel *model, const char *attrname, double *value);
int GRBgetdblattrarray(GRBmodel *model, const char *attrname,
                       int first, int len, double *values);

int GRBaddvar(GRBmodel *model, int numnz, int *vind, double *vval,
              double obj, double lb, double ub, char vtype,
              const char *varname);
int GRBaddconstr(GRBmodel *model, int numnz, int *cind, double *cval,
                 char sense, double rhs, const char *constrname);
int GRBaddqconstr(GRBmodel *model,
                  int numlnz, int *lind, double *lval,
                  int numqnz, int *qrow, int *qcol, double *qval,
                  char sense, double rhs, const char *constrname);

#ifdef __cplusplus
}
#endif
