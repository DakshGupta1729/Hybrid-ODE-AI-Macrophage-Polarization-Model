import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from torch import nn
from torchdiffeq import odeint


STATE_COLUMNS = [
    "V_T",
    "V_T_C1",
    "V_T_T1",
    "V_T_T0",
    "V_T_Mac_M1",
    "V_T_Mac_M2",
    "V_T_IL10",
    "V_T_IL12",
    "V_T_IFNg",
    "V_T_TGFb",
    "V_T_CCL2",
    "V_LN_IL2",
]

IMMUNE_RESIDUAL_COLUMNS = [
    "V_T_T0",
    "V_T_Mac_M1",
    "V_T_Mac_M2",
    "V_T_IL10",
    "V_T_IL12",
    "V_T_IFNg",
    "V_T_TGFb",
    "V_T_CCL2",
]


@dataclass
class Normalizer:
    mean: torch.Tensor
    std: torch.Tensor

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return (x - self.mean) / self.std

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return z * self.std + self.mean


class ResidualDynamics(nn.Module):
    """Neural dx/dt for uncertain Treg/macrophage/cytokine biology."""

    def __init__(self, state_dim: int, param_dim: int, residual_indices: list[int]):
        super().__init__()
        self.residual_indices = residual_indices
        self.net = nn.Sequential(
            nn.Linear(state_dim + param_dim + 1, 96),
            nn.Tanh(),
            nn.Linear(96, 96),
            nn.Tanh(),
            nn.Linear(96, len(residual_indices)),
        )

    def forward(self, t: torch.Tensor, z: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
        squeeze_output = False
        if z.ndim == 1:
            z = z.unsqueeze(0)
            squeeze_output = True
        if p.ndim == 1:
            p = p.unsqueeze(0).expand(z.shape[0], -1)

        t_feature = torch.ones(z.shape[0], 1, device=z.device, dtype=z.dtype) * (t / 400.0)
        residual = self.net(torch.cat([z, p, t_feature], dim=1))
        dz = torch.zeros_like(z)
        dz[:, self.residual_indices] = residual
        if squeeze_output:
            dz = dz.squeeze(0)
        return dz


class PatientODE(nn.Module):
    def __init__(self, dynamics: ResidualDynamics, p: torch.Tensor):
        super().__init__()
        self.dynamics = dynamics
        self.register_buffer("p", p)

    def forward(self, t: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        return self.dynamics(t, z, self.p)


def load_dataset(data_dir: Path, max_patients: int | None):
    trajectories = pd.read_csv(data_dir / "trajectories.csv")
    parameters = pd.read_csv(data_dir / "parameters.csv")

    missing = [col for col in STATE_COLUMNS if col not in trajectories.columns]
    if missing:
        raise ValueError(f"Missing required trajectory columns: {missing}")

    trajectories = trajectories.dropna(subset=STATE_COLUMNS)
    patient_ids = sorted(trajectories["patient_id"].unique())
    if max_patients is not None:
        patient_ids = patient_ids[:max_patients]

    param_columns = [
        col for col in parameters.columns
        if col != "patient_id" and pd.api.types.is_numeric_dtype(parameters[col])
    ]
    if not param_columns:
        raise ValueError("parameters.csv must contain at least one numeric parameter column.")

    patients = []
    for patient_id in patient_ids:
        traj = trajectories[trajectories["patient_id"] == patient_id].sort_values("time")
        param_row = parameters[parameters["patient_id"] == patient_id]
        if traj.empty or param_row.empty:
            continue

        patients.append({
            "patient_id": int(patient_id),
            "time": torch.tensor(traj["time"].to_numpy(dtype=np.float32)),
            "x": torch.tensor(traj[STATE_COLUMNS].to_numpy(dtype=np.float32)),
            "p": torch.tensor(param_row[param_columns].iloc[0].to_numpy(dtype=np.float32)),
        })

    if not patients:
        raise ValueError("No patients with both trajectory and parameter data were found.")

    x_all = torch.cat([patient["x"] for patient in patients], dim=0)
    p_all = torch.stack([patient["p"] for patient in patients], dim=0)
    x_norm = Normalizer(x_all.mean(dim=0), x_all.std(dim=0).clamp_min(1e-6))
    p_norm = Normalizer(p_all.mean(dim=0), p_all.std(dim=0).clamp_min(1e-6))
    return patients, x_norm, p_norm, param_columns


def resolve_data_dir(data_dir: str) -> Path:
    requested = Path(data_dir)
    if (requested / "trajectories.csv").exists() and (requested / "parameters.csv").exists():
        return requested

    kaggle_input = Path("/kaggle/input")
    if kaggle_input.exists():
        matches = [
            path.parent
            for path in kaggle_input.rglob("trajectories.csv")
            if (path.parent / "parameters.csv").exists()
        ]
        if matches:
            print(f"Using Kaggle dataset folder: {matches[0]}")
            return matches[0]

    raise FileNotFoundError(
        "Could not find trajectories.csv and parameters.csv. "
        "Pass --data-dir, or upload a Kaggle dataset containing both files."
    )


def default_output_dir() -> str:
    if Path("/kaggle/working").exists():
        return "/kaggle/working/neural_ode_baseline"
    return "neural_ode/runs/baseline"


def move_normalizer(normalizer: Normalizer, device: torch.device) -> Normalizer:
    return Normalizer(normalizer.mean.to(device), normalizer.std.to(device))


def move_patient(patient: dict, device: torch.device) -> dict:
    return {
        "patient_id": patient["patient_id"],
        "time": patient["time"].to(device),
        "x": patient["x"].to(device),
        "p": patient["p"].to(device),
    }


def train(args):
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device(args.device if args.device != "auto" else ("cuda" if torch.cuda.is_available() else "cpu"))
    print(f"Using device: {device}")
    if device.type == "cuda":
        print(f"GPU: {torch.cuda.get_device_name(0)}")

    data_dir = resolve_data_dir(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    patients, x_norm, p_norm, param_columns = load_dataset(data_dir, args.max_patients)
    patients = [move_patient(patient, device) for patient in patients]
    x_norm = move_normalizer(x_norm, device)
    p_norm = move_normalizer(p_norm, device)

    residual_indices = [STATE_COLUMNS.index(col) for col in IMMUNE_RESIDUAL_COLUMNS]
    dynamics = ResidualDynamics(len(STATE_COLUMNS), len(param_columns), residual_indices).to(device)
    optimizer = torch.optim.AdamW(dynamics.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    print(f"Loaded {len(patients)} patients from {data_dir}")
    print(f"Writing outputs to {out_dir}")

    losses = []
    for epoch in range(1, args.epochs + 1):
        epoch_loss = 0.0
        np.random.shuffle(patients)

        for patient in patients:
            time = patient["time"]
            x = patient["x"]
            p = p_norm.encode(patient["p"])

            z0 = x_norm.encode(x[0])
            target = x_norm.encode(x)
            ode_func = PatientODE(dynamics, p)
            pred = odeint(ode_func, z0, time, method=args.method)

            loss = torch.mean((pred[:, residual_indices] - target[:, residual_indices]) ** 2)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(dynamics.parameters(), 1.0)
            optimizer.step()
            epoch_loss += loss.item()

        epoch_loss /= len(patients)
        losses.append(epoch_loss)
        if epoch == 1 or epoch % args.log_every == 0:
            print(f"epoch={epoch:04d} loss={epoch_loss:.6f}")

    checkpoint = {
        "model_state_dict": dynamics.state_dict(),
        "state_columns": STATE_COLUMNS,
        "residual_columns": IMMUNE_RESIDUAL_COLUMNS,
        "parameter_columns": param_columns,
        "x_mean": x_norm.mean.detach().cpu(),
        "x_std": x_norm.std.detach().cpu(),
        "p_mean": p_norm.mean.detach().cpu(),
        "p_std": p_norm.std.detach().cpu(),
    }
    torch.save(checkpoint, out_dir / "hybrid_neural_ode.pt")

    pd.DataFrame({"epoch": np.arange(1, len(losses) + 1), "loss": losses}).to_csv(
        out_dir / "training_loss.csv", index=False
    )
    plt.figure(figsize=(6, 4))
    plt.plot(losses)
    plt.xlabel("Epoch")
    plt.ylabel("MSE loss")
    plt.tight_layout()
    plt.savefig(out_dir / "training_loss.png", dpi=200)

    plot_example_prediction(patients[0], dynamics, x_norm, p_norm, out_dir, args.method)


def plot_example_prediction(patient, dynamics, x_norm, p_norm, out_dir: Path, method: str):
    dynamics.eval()
    with torch.no_grad():
        time = patient["time"]
        target = patient["x"]
        p = p_norm.encode(patient["p"])
        pred_z = odeint(PatientODE(dynamics, p), x_norm.encode(target[0]), time, method=method)
        pred = x_norm.decode(pred_z)
        time_cpu = time.detach().cpu()
        target_cpu = target.detach().cpu()
        pred_cpu = pred.detach().cpu()

    columns_to_plot = ["V_T_T0", "V_T_Mac_M1", "V_T_Mac_M2", "V_T_IL10", "V_T_IL12", "V_T_TGFb"]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
    for ax, column in zip(axes.ravel(), columns_to_plot):
        idx = STATE_COLUMNS.index(column)
        ax.plot(time_cpu.numpy(), target_cpu[:, idx].numpy(), label="MATLAB", linewidth=2)
        ax.plot(time_cpu.numpy(), pred_cpu[:, idx].numpy(), label="Neural ODE", linestyle="--")
        ax.set_title(column)
        ax.set_xlabel("day")
    axes.ravel()[0].legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"patient_{patient['patient_id']}_fit.png", dpi=200)


def parse_args():
    parser = argparse.ArgumentParser(description="Train a hybrid Neural ODE residual on exported SimBiology trajectories.")
    parser.add_argument("--data-dir", default="neural_ode_data")
    parser.add_argument("--out-dir", default=default_output_dir())
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--max-patients", type=int, default=30)
    parser.add_argument("--method", default="rk4")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--log-every", type=int, default=25)
    return parser.parse_args()


if __name__ == "__main__":
    train(parse_args())
