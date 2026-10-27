# Turn the issued 27-day Kp page into the 3-hour driver file used by the
# NASA Swift challenge. The column layout follows SWORD's ap forecast
# (issue time, bin start, bin end, lead, value) and adds Kp, which is what
# this model actually predicts. Ap is the standard conversion of that Kp.

import argparse
import math
import re
from datetime import datetime, timedelta
from pathlib import Path

# Kp in thirds and the corresponding 3-hour ap (Bartels / GFZ).
KP_AP = [
    (0.0, 0), (1 / 3, 2), (2 / 3, 3), (1.0, 4), (4 / 3, 5), (5 / 3, 6),
    (2.0, 7), (7 / 3, 9), (8 / 3, 12), (3.0, 15), (10 / 3, 18), (11 / 3, 22),
    (4.0, 27), (13 / 3, 32), (14 / 3, 39), (5.0, 48), (16 / 3, 56), (17 / 3, 67),
    (6.0, 80), (19 / 3, 94), (20 / 3, 111), (7.0, 132), (22 / 3, 154), (23 / 3, 179),
    (8.0, 207), (25 / 3, 236), (26 / 3, 300), (9.0, 400),
]


def kp_to_ap(kp):
    if kp <= KP_AP[0][0]:
        return float(KP_AP[0][1])
    if kp >= KP_AP[-1][0]:
        return float(KP_AP[-1][1])
    for (k0, a0), (k1, a1) in zip(KP_AP, KP_AP[1:]):
        if kp <= k1 + 1e-9:
            if a0 == 0:
                return a0 + (kp - k0) / (k1 - k0) * (a1 - a0)
            w = (kp - k0) / (k1 - k0)
            return math.exp(math.log(a0) + w * (math.log(a1) - math.log(a0)))
    return float(KP_AP[-1][1])


def parse_page(html):
    issued = re.search(r"Forecast issued \(UTC\)<b>([^<]+)</b>", html)
    if not issued:
        raise SystemExit("Issue time was not on the page.")
    issue = datetime.strptime(issued.group(1).strip(), "%Y-%m-%d %H:%M:%S")
    start = html.find("<h2>Next 27 days</h2>")
    if start < 0:
        raise SystemExit("The 27-day table was not on the page.")
    table = re.search(r"<table>(.*?)</table>", html[start:], re.S).group(1)
    rows = []
    for raw in re.findall(r"<tr>(.*?)</tr>", table, re.S)[1:]:
        cells = [re.sub("<.*?>", "", c) for c in re.findall(r"<td>(.*?)</td>", raw)]
        rows.append({
            "day": int(cells[0]),
            "end": datetime.strptime(cells[1], "%Y-%m-%d %H:%M:%S"),
            "insitu": float(cells[2]),
            "awsom": float(cells[3]),
            "boost": float(cells[4]),
        })
    return issue, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("html")
    parser.add_argument("out_dir")
    args = parser.parse_args()
    issue, rows = parse_page(Path(args.html).read_text())
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stamp = issue.strftime("%Y%m%dT%H%M")
    three = out / f"UM_Kp_ap_forecast_{stamp}.csv"
    daily_path = out / f"UM_Kp_ap_daily_{stamp}.csv"
    lines = ["issue_time_utc,bin_start_utc,bin_end_utc,lead_h,kp,ap"]
    by_day = {}
    issue_s = issue.strftime("%Y-%m-%d %H:%M")
    for row in rows:
        end = row["end"]
        start = end - timedelta(hours=3)
        lead = int(round((end - issue).total_seconds() / 3600))
        kp = row["boost"]
        ap = kp_to_ap(kp)
        lines.append(
            f"{issue_s},{start.strftime('%Y-%m-%d %H:%M')},{end.strftime('%Y-%m-%d %H:%M')},{lead},{kp:.3f},{ap:.2f}"
        )
        key = start.strftime("%Y-%m-%d")
        by_day.setdefault(key, []).append((kp, ap))
    three.write_text("\n".join(lines) + "\n")
    day_lines = ["date,n_bins,ap_mean,ap_max,kp_mean,kp_max"]
    for day in sorted(by_day):
        vals = by_day[day]
        kps = [v[0] for v in vals]
        aps = [v[1] for v in vals]
        day_lines.append(
            f"{day},{len(vals)},{sum(aps) / len(aps):.2f},{max(aps):.2f},{sum(kps) / len(kps):.3f},{max(kps):.3f}"
        )
    daily_path.write_text("\n".join(day_lines) + "\n")
    same = sum(abs(r["boost"] - r["insitu"]) < 1e-9 for r in rows)
    print(f"issue {issue_s}")
    print(f"bins {len(rows)}  boosted equals in-situ {same}")
    print(f"kp {min(r['boost'] for r in rows):.3f} to {max(r['boost'] for r in rows):.3f}")
    print(three)
    print(daily_path)


if __name__ == "__main__":
    main()
