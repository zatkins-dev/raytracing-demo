from pathlib import Path


def face_id(f: Path):
    return int(f.name.split("_")[1].split(".")[0])


def get_num_faces(data_dir: Path, base_name: str) -> int:
    return max(map(face_id, data_dir.glob(f"{base_name}_*.plt")))


if __name__ == '__main__':
    print(get_num_faces(Path(__file__).parent.parent / "data", "hit"))
