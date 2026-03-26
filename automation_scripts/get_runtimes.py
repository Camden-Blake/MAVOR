import csv
import re
from pathlib import Path

logfile = Path("log.log")
outfile = Path("summary.csv")

material_re = re.compile(r"Processing\s+(tsl-\S+)\.")
step_re = re.compile(r"--Completed in ([0-9.]+) seconds")
total_re = re.compile(r"Time to generate OTF file for (\S+):\s+([0-9.]+) seconds")

STEP_ORDER = [
    "LEAPR",
    "Linearize SAB",
    "Union Energy",
    "Unioned SAB",
    "OTF",
    "Elastic",
]

def parse_log(text):
    rows = []
    current_mat = None
    step_times = []
    lines = text.splitlines()
    for line in lines:
        mat_match = material_re.search(line)
        if mat_match:
            current_mat = mat_match.group(1)
            step_times = []
            continue
        step_match = step_re.search(line)
        if step_match:
            step_times.append(round(float(step_match.group(1)), 2))
            continue
        total_match = total_re.search(line)
        if total_match:
            mat_name = total_match.group(1)
            total_time = round(float(total_match.group(2)), 2)
            if len(step_times) != len(STEP_ORDER):
                raise RuntimeError(
                    f"Material {mat_name} has {len(step_times)} steps, "
                    f"expected {len(STEP_ORDER)}."
                )
            row = {
                "Material": mat_name,
                **{step: t for step, t in zip(STEP_ORDER, step_times)},
                "Total": total_time
            }
            rows.append(row)
    return rows


def write_csv(rows, outfile):
    if not rows:
        print("No data parsed.")
        return
    fieldnames = ["Material"] + STEP_ORDER + ["Total"]
    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    text = logfile.read_text()
    rows = parse_log(text)
    write_csv(rows, outfile)
    print(f"Wrote {outfile} with {len(rows)} materials.")