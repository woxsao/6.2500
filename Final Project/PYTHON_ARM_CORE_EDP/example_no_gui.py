from edp import ARMCore

# This is an example of how to use the ARMCore simulator without the
# jupyter notebook GUI.

# Create a re-usable instance of the ARMCore simulator. This will load
# the ARM core circuit model and layout into memory and do some
# pre-calculations.
armcore = ARMCore()

# we can re-use this instance to calculate the delay and energy
# simple calculation:
results = armcore.calculate_energy_delay(
    i_on=600, # in uA/um
    i_off=400,   # in nA/um
    v_dd=1.8,  # supply voltage
    c_gs=1.0,     # normalized to 1
)

print("RESULTS:")
for key, value in results.items():
    print(f"{key} = {value}")


# we can loop values and re-use the instance:
# for example, here we can see effect of reducing capacitance c_gs
for (i_on, i_off, v_dd, c_gs) in [
    (600, 400, 1.8, 1.0),
    (600, 400, 1.8, 0.5),
    (600, 400, 1.8, 0.1),
]:
    # calculate edp
    results = armcore.calculate_energy_delay(
        i_on=i_on,
        i_off=i_off,
        v_dd=v_dd,
        c_gs=c_gs,
    )

    print("====================================================================")
    print(f"i_on={i_on}, i_off={i_off}, v_dd={v_dd}, c_gs={c_gs}")
    print("RESULTS:")
    for key, value in results.items():
        print(f"{key} = {value}")
