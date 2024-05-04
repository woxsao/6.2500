; --- Define parameters and list values ---

; filename
(define filename "pn")

; device dimensions (units in um)
; p-type length lp and ntype side ln
(define LP 1)
(define LN 1)
; thickness
(define LY 0.5)

; doping concentrations (cm^-3)
(define P_CONC 1E16)
(define N_CONC 1E16)

; meshing refinement size, sentaurus does finite element analysis
; this sizing determines how fine the grid is. finer grids take
; longer to simulate but can be more spatially accurate
; (units in um)
(define XMESH_MIN 0.005)
(define XMESH_MAX 0.05)
(define YMESH_MIN 0.01)
(define YMESH_MAX 0.05)

; --- derived variables ---
(define xp (- 0 LP))     ; p-contact edge xp, makes it negative
(define xn LN)           ; n-contact edge xn
(define ytop LY)         ; y value "top" of device
(define ymid (* 0.5 LY)) ; y value middle of device

; --- Program run ---

(begin 
(sde:clear) 

; --- Create structure ---
; using a single rectangle of silicon, we will specify doping regions later
(sdegeo:create-rectangle (position xp 0 0) (position xn ytop 0) "Silicon" "Bulk")

; --- Contacts ---

; Contact definitions and assignments

; P-contact
(sdegeo:define-contact-set "anode" 4.0  (color:rgb 1.0 0.0 0.0 ) "##")
(sdegeo:set-current-contact-set "anode") 
(sdegeo:define-2d-contact (find-edge-id (position xp ymid 0)) "anode")

; N-contact
(sdegeo:define-contact-set "cathode" 4.0  (color:rgb 0.0 0.0 1.0 ) "==")
(sdegeo:set-current-contact-set "cathode")
(sdegeo:define-2d-contact (find-edge-id (position xn ymid 0)) "cathode")

; --- Reference windows for doping and meshing ---

(sdedr:define-refeval-window "RW_Bulk" "Rectangle" (position xp 0 0) (position xn ytop 0))
(sdedr:define-refeval-window "RW_P" "Rectangle" (position xp 0 0) (position 0 ytop 0))
(sdedr:define-refeval-window "RW_N" "Rectangle" (position 0 0 0) (position xn ytop 0))

; Doping

(sdedr:define-constant-profile "C_P" "BoronActiveConcentration" P_CONC)
(sdedr:define-constant-profile "C_N" "ArsenicActiveConcentration" N_CONC)

(sdedr:define-constant-profile-placement "C_PL_P" "C_P" "RW_P" )
(sdedr:define-constant-profile-placement "C_PL_N" "C_N" "RW_N" )

; --- Creating Mesh ---

(sdedr:define-refinement-size "M_Bulk" XMESH_MAX YMESH_MAX XMESH_MIN YMESH_MIN)

(sdedr:define-refinement-placement "M_PL_Bulk" "M_Bulk" "RW_Bulk" )

; save model
(sde:save-model "si_pn_diode")

; Build Mesh
(sde:build-mesh "snmesh" "-d" filename)
)

