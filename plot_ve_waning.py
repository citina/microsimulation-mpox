import numpy as np
import matplotlib.pyplot as plt

# Constants
vac1_normal = 0.721
vac2_normal = 0.878
vac1_plwh = 0.51
vac2_plwh = 0.702

# Time points (weeks)
weeks = np.arange(0, 65)  # 0 to 64 weeks

# Previous version (waning to 50%)
def previous_ve_normal(week):
    if week <= 2:
        return 0.3605 * week
    elif week <= 4:
        return 0.721
    elif week <= 6:
        return max(round(-0.014 * week + 0.777, 3), 0)
    elif week <= 58:
        return max(round(-0.017 * week + 0.986, 3), 0)
    else:
        return 0.439

def previous_ve_plwh(week):
    if week <= 2:
        return 0.255 * week
    elif week <= 4:
        return 0.51
    elif week <= 6:
        return max(round(-0.01 * week + 0.56, 3), 0)
    elif week <= 58:
        return max(round(-0.014 * week + 0.786, 3), 0)
    else:
        return 0.351

# Current version (waning to 50%)
def current_ve_normal(week):
    if week <= 2:
        return 0.3605 * week
    elif week <= 4:
        return 0.721
    elif week <= 6:
        return max(round(-0.014 * week + 0.777, 3), 0)
    elif week <= 58:
        return max(round(-0.00844 * (week - 6) + 0.878, 3), 0)
    else:
        return 0.439

def current_ve_plwh(week):
    if week <= 2:
        return 0.255 * week
    elif week <= 4:
        return 0.51
    elif week <= 6:
        return max(round(-0.01 * week + 0.56, 3), 0)
    elif week <= 58:
        return max(round(-0.00675 * (week - 6) + 0.702, 3), 0)
    else:
        return 0.351

# Calculate VE values
prev_normal = [previous_ve_normal(w) for w in weeks]
prev_plwh = [previous_ve_plwh(w) for w in weeks]
curr_normal = [current_ve_normal(w) for w in weeks]
curr_plwh = [current_ve_plwh(w) for w in weeks]

# Create the plot
plt.figure(figsize=(12, 6))

# Plot previous version
plt.plot(weeks, prev_normal, 'b-', label='Previous - Normal')
plt.plot(weeks, prev_plwh, 'b--', label='Previous - PLWH')

# Plot current version
plt.plot(weeks, curr_normal, 'r-', label='Current - Normal')
plt.plot(weeks, curr_plwh, 'r--', label='Current - PLWH')

# Add labels and title
plt.xlabel('Weeks since first dose')
plt.ylabel('Vaccine Effectiveness')
plt.title('VE Waning Over Time (50% waning case)')
plt.grid(True)
plt.legend()

# Add vertical lines for key time points
plt.axvline(x=2, color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=4, color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=6, color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=58, color='gray', linestyle='--', alpha=0.5)

# Add text labels for key time points
plt.text(2, 0.1, 'Week 2', rotation=90, alpha=0.5)
plt.text(4, 0.1, 'Week 4', rotation=90, alpha=0.5)
plt.text(6, 0.1, 'Week 6', rotation=90, alpha=0.5)
plt.text(58, 0.1, 'Week 58', rotation=90, alpha=0.5)

# Save the plot
plt.savefig('ve_waning_comparison.png')
plt.close() 