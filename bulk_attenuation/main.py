from dataclasses import dataclass
import unyt as u


@dataclass(frozen=True)
class PhaseProperties:
    """Holds the physical and thermodynamic properties of a single phase."""

    name: str
    bulk_modulus: u.unyt_quantity  # Dimensions: pressure
    thermal_expansion: u.unyt_quantity  # Dimensions: 1/temperature
    specific_heat: u.unyt_quantity  # Dimensions: energy / (mass * temperature)
    density: u.unyt_quantity  # Dimensions: mass / volume


def calculate_adiabatic_temperature_gradient(
    phase: PhaseProperties, T_0: u.unyt_quantity
) -> u.unyt_quantity:
    """Calculates the adiabatic temperature change per unit pressure change (dT/dP).

    Formula: dT/dP = (alpha * T0) / (rho * Cp)
    """
    numerator = phase.thermal_expansion * T_0
    denominator = phase.density * phase.specific_heat
    return numerator / denominator


def calculate_temperature_change(
    phase: PhaseProperties,
    T_0: u.unyt_quantity,
    dP: u.unyt_quanity, # dimensions: pressure
) -> u.unyt_quantity:
    return calculate_adiabatic_temperature_gradient(phase, T_0) * dP



# --- 1. Define the Environmental Baseline ---
# Ambient Upper Mantle Temperature
T0 = 1600.0 * u.K

# --- 2. Define the Phase Endmembers ---
# Solid Matrix (Olivine-dominated)
solid_matrix = PhaseProperties(
    name="Solid Matrix (Olivine)",
    bulk_modulus=130.0 * u.GPa,
    thermal_expansion=3e-5 / u.K,
    specific_heat=1200.0 * u.J / (u.kg * u.K),
    density=3300.0 * u.kg / (u.m**3),
)

# Fluid Phase (Basaltic Melt)
basaltic_melt = PhaseProperties(
    name="Fluid Phase (Basaltic Melt)",
    bulk_modulus=15.0 * u.GPa,
    thermal_expansion=6e-5 / u.K,
    specific_heat=1500.0 * u.J / (u.kg * u.K),
    density=2700.0 * u.kg / (u.m**3),
)

# --- 3. Run Calculations ---
dT_dP_solid = calculate_adiabatic_temperature_gradient(solid_matrix, T0)
dT_dP_fluid = calculate_adiabatic_temperature_gradient(basaltic_melt, T0)

# Quantify the absolute gradient mismatch across the boundary
dT_dP_diff = abs(dT_dP_fluid - dT_dP_solid)

# --- 4. Display Results ---
print(f"Mantle Baseline Temperature: {T0}\n")

print(f"{solid_matrix.name}:")
print(f"  dT/dP = {dT_dP_solid.to('K/Pa'):.3e}")
print(f"  dT/dP = {dT_dP_solid.to('K/MPa'):.4f}\n")

print(f"{basaltic_melt.name}:")
print(f"  dT/dP = {dT_dP_fluid.to('K/Pa'):.3e}")
print(f"  dT/dP = {dT_dP_fluid.to('K/MPa'):.4f}\n")

print("-" * 50)
print(f"Induced Interphase Thermal Gradient Mismatch:")
print(f"  Δ(dT/dP) = {dT_dP_diff.to('K/Pa'):.3e}")
print(f"  Δ(dT/dP) = {dT_dP_diff.to('K/MPa'):.4f}")

def main():
    print("Hello from bulk-attenuation!")


if __name__ == "__main__":
    main()
