;; Setting Parameters:


;; Defined Parameters:


;; Contact Sets:

;; Work Planes:
(sde:workplanes-init-scm-binding)

;; Defined ACIS Refinements:
(sde:refinement-init-scm-binding)

;; Restore GUI session parameters:
(sde:set-window-position 179 82)
(sde:set-swindow-size 840 800)
(sde:set-window-style "Windows")
(sde:set-background-color 0 127 178 204 204 204)
(sde:scmwin-set-prefs "DejaVu Sans Mono" "Normal" 10 100 )

;; NDK: define structure here!
(sde:clear)

(sdegeo:set-default-boolean "ABA")

;; assign geometries
;; silicon channel:
(sdegeo:create-rectangle (position 0 0 0.0 )  (position 0.035 0.035 0.0 ) "Silicon" "channel" )
;; gate oxide:
(sdegeo:create-rectangle (position 0 0 0 )  (position 0.035 -0.004 0 ) "SiO2" "gate_oxide" )
;; Gate metal electrode:
(sdegeo:create-rectangle (position 0 -0.004 0 )  (position 0.035 -0.05 0 ) "Aluminum" "gate_electrode" )
;; spacer oxide between metal gate and metal source:
(sdegeo:create-rectangle (position 0 0 0 )  (position -0.035 -0.05 0 ) "SiO2" "spacerL" )
;; spacer oxide between metal gate and metal drain:
(sdegeo:create-rectangle (position 0.035 0 0 )  (position 0.07 -0.05 0 ) "HfO2" "spacerR" )
;; Drain metal electrode:
(sdegeo:create-rectangle (position 0.07 0 0 )  (position 0.105 -0.05 0.0 ) "Aluminum" "drain_electrode" )
;; Source metal electrode:
(sdegeo:create-rectangle (position -0.07 0 0 )  (position -0.035 -0.05 0 ) "Aluminum" "source_electrode" )
;; Body metal electrode:
(sdegeo:create-rectangle (position -0.07 0.6 0 )  (position 0.105 0.5 0 ) "Aluminum" "body_electrode" )
;; silicon under source metal:
(sdegeo:create-rectangle (position -.07 .02 0 )  (position -0.035 0 0 ) "Silicon" "source_n" )
;; silicon between silicon under source metal and silicon channel, called source extension:
(sdegeo:create-rectangle (position -.035 .01 0 )  (position 0 0 0 ) "Silicon" "source_n_ext" )
;; silicon between silicon under drain metal and silicon channel, called drain extension:
(sdegeo:create-rectangle (position .035 .01 0 )  (position .07 0 0 ) "Silicon" "drain_n_ext" )
;; silicon under drain metal:
(sdegeo:create-rectangle (position .07 .02 0 )  (position .105 0 0 ) "Silicon" "drain_n" )
;; body silicon:
(sdegeo:create-rectangle (position -.07 .5 0 )  (position .105 .035 0 ) "Silicon" "body" )



;; assign doping
;; channel doping:
(sdedr:define-constant-profile "constant_channel_doping" "BoronActiveConcentration" 3e18)
(sdedr:define-constant-profile-region "constant_channel_doping_placement" "constant_channel_doping" "channel")
;; drain extension doping:
(sdedr:define-constant-profile "constant_drain_ext_doping" "PhosphorusActiveConcentration" 5e+18)
(sdedr:define-constant-profile-region "constant_drain_ext_doping_placement" "constant_drain_ext_doping" "drain_n_ext")
;; drain doping:
(sdedr:define-constant-profile "constant_drain_doping" "PhosphorusActiveConcentration" 5e+18)
(sdedr:define-constant-profile-region "constant_drain_doping_placement" "constant_drain_doping" "drain_n")
;; source extension doping:
(sdedr:define-constant-profile "constant_source_ext_doping" "PhosphorusActiveConcentration" 5e+19)
(sdedr:define-constant-profile-region "constant_source_ext_doping_placement" "constant_source_ext_doping" "source_n_ext")
;; source doping:
(sdedr:define-constant-profile "constant_source_doping" "PhosphorusActiveConcentration" 5e+19)
(sdedr:define-constant-profile-region "constant_source_doping_placement" "constant_source_doping" "source_n")
;; body doping:
 (sdedr:define-constant-profile "constant_body_doping" "BoronActiveConcentration" 5e+19)
(sdedr:define-constant-profile-region "constant_body_doping_placement" "constant_body_doping" "body")



;; define contacts
(sdegeo:define-contact-set "Gate_contact" 4  (color:rgb 1 0 0 ) "##" )
(sdegeo:define-contact-set "Body_contact" 4  (color:rgb 1 1 0 ) "##" )
(sdegeo:define-contact-set "Source_contact" 4  (color:rgb 1 0 1 ) "##" )
(sdegeo:define-contact-set "Drain_contact" 4  (color:rgb 1 1 1 ) "##" )

;; assign location of contacts
(sdegeo:set-current-contact-set "Gate_contact")
(sdegeo:set-contact-boundary-edges (list (car (find-body-id (position 0.02 -0.04 0)))) "Gate_contact")
(sdegeo:set-current-contact-set "Source_contact")
(sdegeo:set-contact-boundary-edges (list (car (find-body-id (position -.06 -0.04 0)))) "Source_contact")
(sdegeo:set-current-contact-set "Drain_contact")
(sdegeo:set-contact-boundary-edges (list (car (find-body-id (position 0.1 -0.04 0)))) "Drain_contact")
(sdegeo:set-current-contact-set "Body_contact")
(sdegeo:set-contact-boundary-edges (list (car (find-body-id (position 0.02 0.55 0)))) "Body_contact")


;; make mesh
;; define and set top mesh, everything above the surface of the silicon:
(sdedr:define-refeval-window "RefWin_1" "Rectangle"  (position -0.07 0 0) (position 0.105 -0.05 0))
(sdedr:define-refinement-size "RefinementDefinition_1" 0.002 0.001 0.005 0.002 )
(sdedr:define-refinement-placement "RefinementPlacement_1" "RefinementDefinition_1" "RefWin_1" )
;; define and set the bottom mesh, everything below the surface of the silicon:
(sdedr:define-refeval-window "RefWin_2" "Rectangle"  (position -0.07 .6 0) (position 0.105 0 0))
(sdedr:define-multibox-size "MultiboxDefinition_1" 0.05 0.03 .03 .0002 1 1.35 )
(sdedr:define-multibox-placement "MultiboxPlacement_1" "MultiboxDefinition_1" "RefWin_2" )
;; define and set the mesh over the channel:
(sdedr:define-refinement-size "RefinementDefinition_1" 0.002 .005 .001 0.002 )
(sdedr:define-refinement-region "RefinementPlacement_3" "RefinementDefinition_1" "channel" )


;; save file here
;; MANUALLY SAVE THE FILE WITH A FILENAME HERE BEFORE DOING MESHING!!!

(sde:save-model "Si_MOSFET")

;; manually build mesh
;; mesh->build mesh

;; save file

(sde:build-mesh "snmesh" "" "Mesh_Si_MOSFET")
