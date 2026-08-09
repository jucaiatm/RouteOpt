#include "gurobi_c.h"

struct GRBenv {};
struct GRBmodel {};

extern "C" {

const char *GRBgeterrormsg(GRBenv *) {
    return "Gurobi stub: solver routines are unavailable in this test";
}

GRBenv *GRBgetenv(GRBmodel *) {
    static GRBenv env;
    return &env;
}

int GRBloadenv(GRBenv **env, const char *) {
    if (env != nullptr) *env = new GRBenv;
    return 0;
}

void GRBfreeenv(GRBenv *env) {
    delete env;
}

int GRBnewmodel(GRBenv *, GRBmodel **model, const char *, int,
                double *, double *, double *, char *, char **) {
    if (model != nullptr) *model = new GRBmodel;
    return 0;
}

void GRBfreemodel(GRBmodel *model) {
    delete model;
}

int GRBupdatemodel(GRBmodel *) { return 0; }
int GRBoptimize(GRBmodel *) { return 0; }
int GRBsetintparam(GRBenv *, const char *, int) { return 0; }
int GRBsetintattr(GRBmodel *, const char *, int) { return 0; }
int GRBsetdblattrelement(GRBmodel *, const char *, int, double) { return 0; }

int GRBgetintattr(GRBmodel *, const char *, int *value) {
    if (value != nullptr) *value = 0;
    return 0;
}

int GRBgetdblattr(GRBmodel *, const char *, double *value) {
    if (value != nullptr) *value = 0.0;
    return 0;
}

int GRBgetdblattrarray(GRBmodel *, const char *, int, int len, double *values) {
    if (values != nullptr) {
        for (int i = 0; i < len; ++i) values[i] = 0.0;
    }
    return 0;
}

int GRBaddvar(GRBmodel *, int, int *, double *, double, double, double, char, const char *) {
    return 0;
}

int GRBaddconstr(GRBmodel *, int, int *, double *, char, double, const char *) {
    return 0;
}

int GRBaddqconstr(GRBmodel *, int, int *, double *, int, int *, int *, double *,
                  char, double, const char *) {
    return 0;
}

} // extern "C"
