import os
import csv
import matplotlib.pyplot as plt
import numpy as np

data_dir = "./results"
otimizado = "v2"
original = "v1"
os.makedirs("graphs", exist_ok=True)

allData = {
    "gradient": {
        "otimizado": {},
        "original": {}
    },
    "residue": {
        "otimizado": {},
        "original": {}
    }
}

n_vals_str = ["32", "64", "128", "256", "512", "1000"] # , "2000", "4000", "8000", "9000", "10000", "20000"]

n_vals = np.array(["32", "64", "128", "256", "512", "1000" ]) #, "2000", "4000", "8000", "9000", "10000", "20000"])

n_vals_int = np.array([int(n) for n in n_vals_str])

OP_MAP = {
    "gradient": {"id": "op1", "label": "Gradiente", "color_opt": "blue", "color_nonopt": "orange"},
    "residue": {"id": "op2", "label": "Resíduo", "color_opt": "green", "color_nonopt": "red"}
}

metrics = []
for file in os.listdir(f"{data_dir}/{otimizado}"):
    filename = os.fsdecode(file)
    if filename.endswith(".csv"):
        metrics.append(filename)

        for op in allData.keys():
            for ver in allData[op].keys():
                allData[op][ver][filename] = {}

        with open(f"{data_dir}/{otimizado}/{filename}", mode ='r') as f:
            next(f)
            csvFile = csv.reader(f)
            for lines in csvFile:
                N = lines[0]
                allData["gradient"]["otimizado"][filename][N] = float(lines[1])
                allData["residue"]["otimizado"][filename][N] = float(lines[2])

        with open(f"{data_dir}/{original}/{filename}", mode ='r') as f:
            next(f)
            csvFile = csv.reader(f)
            for lines in csvFile:
                N = lines[0]
                allData["gradient"]["original"][filename][N] = float(lines[1])
                allData["residue"]["original"][filename][N] = float(lines[2])

for operation, op_info in OP_MAP.items():
    div_factor = 1

    if operation == "gradient":
        div_factor = 25

    for metric in metrics:
        metric_title = str.title(metric[0:-4].replace("_", " "))

        y_label = metric_title

        if metric_title == "Exec Time":
            y_label += " (ms)"
        elif metric_title == "L2Cache":
            y_label += " Miss Ratio"
        elif metric_title == "L3":
            y_label += " MBytes/s"
        else:
            y_label += "/s"

        raw_values_opt    = [allData[operation]["otimizado"]   [metric][n] / div_factor for n in n_vals_str]
        raw_values_nonopt = [allData[operation]["original"][metric][n] / div_factor for n in n_vals_str]

        fig, axs = plt.subplots(1, 1, figsize=(12, 6))

        fig.suptitle(f"{metric_title} por Tamanho do SL - {op_info['label']} ", fontsize=14)

        axs.plot(n_vals_int, raw_values_opt,
                 marker='o', linestyle='-', color=op_info['color_opt'],
                 label ='v2 - Otimizada')

        axs.plot(n_vals_int, raw_values_nonopt,
                 marker='x', linestyle='--', color=op_info['color_nonopt'],
                 label ='v1 - Não Otimizada')

        axs.set_ylabel(y_label, fontsize=12)

        if metric_title == "Exec Time":
            axs.set_yscale('log')
            axs.set_title("Eixo Y em Escala Logarítmica para Visualizar Diversas Ordens de Magnitude", fontsize=10, color='red')

        axs.set_xlabel('Tamanho do Sistema Linear (N)', fontsize=12)
        axs.set_xscale('log')
        axs.set_xticks(n_vals_int)
        axs.set_xticklabels(n_vals, fontsize=8)

        # axs.grid(True, which="both", ls="--", linewidth=0.5)
        axs.legend(title="Versão do Código")

        plt.tight_layout(rect=(0.0, 0.03, 1.0, 0.95))

        filename = f"{op_info['id']}_{operation}_{metric_title.replace(' ', '_')}.png"
        plt.savefig(f"graphs/{filename}", dpi=300)

        plt.close(fig)
