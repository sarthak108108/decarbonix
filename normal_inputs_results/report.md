# CarbonOS results

Coal and biomass continuously ON; minimum-on requirements satisfied; starts=0; 1h intervals.

Gate flags compare to the official normal thresholds; they are not required for the disturbance case.

```json
{
  "cost": 5045753.333333184,
  "emissions": 856.307333333292,
  "cost_saving_pct": 8.702569818388051,
  "carbon_saving_pct": 13.240701714573966,
  "coal_t": 297.5733333333159,
  "biomass_t": 220.0,
  "gas_sm3": 0.0,
  "steam_totals": [
    1653.1851851850881,
    814.8148148148148,
    0.0,
    329.0
  ],
  "solar_mwh": 74.0,
  "grid_mwh": 201.70000000000005,
  "flex_mwh": 18.0,
  "steam_emissions": 721.297333333292,
  "grid_emissions": 135.01,
  "steam_cost": 3739373.333333184,
  "grid_cost": 1306380.0,
  "carbon_gate_pass": true,
  "cost_gate_pass": true
}
```

## Feasibility checks

```json
{
  "steam_min": -9.740119821799453e-11,
  "steam_max": 2.519999999999996,
  "electric_balance_error": -0.0,
  "coal_bounds": 13.415343915350888,
  "bio_bounds": 0.0,
  "gas_bounds": 0.0,
  "waste_bounds": 0.0,
  "solar_bounds": 0.0,
  "grid_bounds": 5.5,
  "flex_bounds": 0.0,
  "flex_forbidden": -0.0,
  "flex_daily_error": -0.0,
  "continuity": 0.0,
  "coal_fuel_margin": 52.42666666668413,
  "bio_fuel_margin": 0.0,
  "gas_fuel_margin": 20000.0,
  "coal_ramp_margin": 1.0000000000000142,
  "bio_ramp_margin": -1.4210854715202004e-14,
  "gas_ramp_margin": 45.0,
  "waste_ramp_margin": 17.0
}
```

Full hourly set-points, costs, emissions and binding constraints: hourly_schedule.csv.
