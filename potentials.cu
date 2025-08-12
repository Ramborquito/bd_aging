#include "encabezados.h"

//Note: be careful with signs.

// =============================================== 2023-05-26 ========================
// HOST
// ===================================================================================

void potential_wca49_50_hst(float sigma, float dist_inv, float invT,float AA,
                            float &potential, float &normal_force){
    
    float sor = dist_inv*sigma;
    float sor6;
    sor6 = sor*sor*sor;
    sor6 = sor6*sor6;
    float sor12=sor6*sor6;
    float sor24=sor12*sor12;
    float sor48=sor24*sor24;
    float sor50=sor48*sor*sor;
    
    potential += EPS*AA*invT*(sor50 - sor48*sor)+EPS*invT;
    normal_force = -EPS*AA*invT*(50*sor50-49*sor48*sor)*dist_inv;
}

// ===================================================================================

void potential_wca49_50_AO_hst(float sigma, float dist_inv, float invT,float AA, float sigma_pol, float phi_pol,
                               float &potential, float &normal_force){
    float size_ratio=sigma/sigma_pol;
    float sor = dist_inv*sigma;
    float ros = 1.0/sor;
    
    
    if(ros<=1.0){
        float sor6 = sor*sor*sor*sor*sor*sor;
        float sor12=sor6*sor6;
        float sor24=sor12*sor12;
        float sor48=sor24*sor24;
        float sor50=sor48*sor*sor;
        
        potential += EPS*AA*invT*(sor50 - sor48*sor)+EPS*invT;
        normal_force = -EPS*AA*invT*(50*sor50-49*sor48*sor)*dist_inv;
    }
    else{
        float ros3 = ros*ros*ros;
        float eta_1 = size_ratio+1.0;
        float eta2 = eta_1*eta_1;
        potential += -EPS*phi_pol*(eta2*eta_1-1.5*eta2*ros+0.5*ros3);
        normal_force = EPS*phi_pol*(1.5*eta2*ros-1.5*ros3)*dist_inv;
    }
    
}

void potential_wca_hst(float sigma, float dist_inv, float &potential, float &normal_force){
    float sor = dist_inv*sigma;
    float sor6, sor12;
    sor6 = sor*sor*sor;
    sor6 = sor6*sor6;
    sor12=sor6*sor6;

    potential += EPS * 4 * (sor12 - sor6) + 1;
    normal_force = - EPS * 4 * (12 * sor12 - 6 * sor6) * dist_inv;

}
// ===================================================================================
// DEVICE
// ===================================================================================

__device__ void potential_wca49_50_dev(float sigma, float dist_inv, float invT,float AA,
                                       float &potential, float &normal_force){
    
    float sor = dist_inv*sigma;
    
    float sor6 = sor*sor*sor*sor*sor*sor;
    float sor12=sor6*sor6;
    float sor24=sor12*sor12;
    float sor48=sor24*sor24;
    float sor50=sor48*sor*sor;
    
    potential += EPS*AA*invT*(sor50 - sor48*sor)+EPS*invT;
    normal_force = -EPS*AA*invT*(50*sor50-49*sor48*sor)*dist_inv;
}

__device__ void potential_wca49_50_AO_dev(float sigma, float dist_inv, float invT,float AA, float sigma_pol, float phi_pol,
                                          float &potential, float &normal_force){
    float size_ratio=sigma/sigma_pol;
    float sor = dist_inv*sigma;
    float ros = 1.0/sor;
    
    
    if(ros<=1.0){
        float sor6 = sor*sor*sor*sor*sor*sor;
        float sor12=sor6*sor6;
        float sor24=sor12*sor12;
        float sor48=sor24*sor24;
        float sor50=sor48*sor*sor;
        
        potential += EPS*AA*invT*(sor50 - sor48*sor)+EPS*invT;
        normal_force = -EPS*AA*invT*(50*sor50-49*sor48*sor)*dist_inv;
    }
    else{
        float ros3 = ros*ros*ros;
        float eta_1 = size_ratio+1.0;
        float eta2 = eta_1*eta_1;
        potential += -EPS*phi_pol*(eta2*eta_1-1.5*eta2*ros+0.5*ros3);
        normal_force = EPS*phi_pol*(1.5*eta2*ros-1.5*ros3)*dist_inv;
    }
    
}


__device__ void potential_wca_dev(float sigma, float dist_inv, float &potential, float &normal_force){
    float sor = dist_inv*sigma;
    float sor6, sor12;
    sor6 = sor*sor*sor;
    sor6 = sor6*sor6;
    sor12=sor6*sor6;

    potential += EPS * 4 * (sor12 - sor6) + 1;
    normal_force = - EPS * 4 * (12 * sor12 - 6 * sor6) * dist_inv;

}