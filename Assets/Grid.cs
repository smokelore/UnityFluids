using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Cell {
	public int[] index = new int[3];
	public Vector3 position;
	public Vector3 pressure;
	public Vector3 velocity;
	public Vector3 prev_velocity;
	public float[] density = new float[3];
	public float[] prev_density = new float[3];
	public bool isBoundary = false;
	public ParticleSystem.Particle particle;

	public Vector3 allowVelocity = Vector3.one;

	public void SetVelocity(Vector3 new_velocity) {
		this.prev_velocity = this.velocity;
		//new_velocity = Vector3.ClampMagnitude(new_velocity, 1f);
		this.velocity = new_velocity;
		//Debug.DrawLine(position, position + velocity, Color.red);
	}

	public void ChangeVelocity(Vector3 delta) {
		SetVelocity(1f * this.velocity + 1f * delta);
	}

	public void SetAsBoundary() {
		this.isBoundary = true;
		//Debug.Log("" + index[0] + " " + index[1] + " " + index[2] + " : " + isBoundary);
	}

	public void SetDensity(float new_density, int colorIndex) {
		float[] new_density_rgb = this.density;
		new_density_rgb[colorIndex] = new_density;

		this.prev_density = this.density;
		this.density = new_density_rgb;

		UpdateColor();
	}

	public void ChangeDensity(float delta_density, int colorIndex) {
		SetDensity(this.density[colorIndex] + delta_density, colorIndex);
	}

	public void UpdateColor() {
		if (!isBoundary) {
			Color rgb = new Color(density[0], density[1], density[2], 1f);
			if (!Grid.instance.UseTexture) {
				rgb.a = Mathf.Max(density[0], density[1], density[2]);
			}
			particle.color = rgb;
		}
	}

	public bool PointIsWithin(Vector3 testPosition) {
		if (!isBoundary) {
			Vector3 minExtents = this.position - Vector3.one * this.particle.size/2;
			Vector3 maxExtents = this.position + Vector3.one * this.particle.size/2;
			return (testPosition.x < maxExtents.x && testPosition.x >= minExtents.x
					&& testPosition.y < maxExtents.y && testPosition.y >= minExtents.y
					&& testPosition.z < maxExtents.z && testPosition.z >= minExtents.z);
		} else {
			return false;
		}
	}
}






public class Source {
	public int[] index = new int[3];
	public Vector3 direction;
	public float amount;

	public float initAmount;
	private float initTime;
	private float duration;

	private int colorIndex;

	public Source(int[] index, Vector3 direction, float amount, int colorIndex) {
		Initialize(index, direction, amount, colorIndex);
	}

	private void Initialize(int[] index, Vector3 direction, float amount, int colorIndex) {
		this.index = index;
		this.direction = Vector3.ClampMagnitude(direction, 1f);
		this.colorIndex = colorIndex;
		this.initAmount = amount;
		this.initTime = Time.time;
		this.duration = Grid.instance.SourceDuration;

		/*Grid.instance.cells[index[0], index[1], index[2]].ChangeDensity(amount * 10 * Grid.instance.dt);
		Grid.instance.cells[index[0], index[1], index[2]].ChangeVelocity(amount * direction * Grid.instance.dt);
		Grid.instance.UpdateParticle(index[0], index[1], index[2]);*/
		//Grid.instance.cells[index[0], index[1], index[2]].ChangeVelocity(amount * direction);
	}

	public void UpdateEmission() {
		float progress = (Time.time - initTime)/duration;
		this.amount = Mathf.Lerp(0f, initAmount, Mathf.Pow(progress+0.25f, 2f));

		if (!Grid.instance.UseTexture) {
			Grid.instance.cells[index[0], index[1], index[2]].ChangeDensity(amount * 10 * Grid.instance.dt, colorIndex);
		}
		Grid.instance.cells[index[0], index[1], index[2]].ChangeVelocity(amount * direction * Grid.instance.dt);
		Grid.instance.UpdateParticle(index[0], index[1], index[2]);
	}

	public bool isAlive() {
		return (Time.time - initTime < duration);
	}
}
















public class Grid : MonoBehaviour {
	public static Grid instance;
	public ParticleSystem system;
	public Vector3 GRID_RESOLUTION; 
	public Cell[,,] cells;

	private Vector3 ORIGIN;
	public float GRID_SIZE;
	private int[] RESOLUTION = new int[3];
	private float PARTICLE_SIZE;

	public float dt = 0;
	public float DIFFUSION_RATE;

	public List<Source> sources = new List<Source>();
	private ParticleSystem.Particle[] particles;

	public GameObject UserObject;
	public float SourceStrength;
	public float SourceDuration;

	public int relaxationIterations;

	public int ChosenColor;

	public bool UseTexture;
	public Texture2D ImportedTexture;

	public bool ExternalForce;
	public Vector3 externalForce;

	void Start() {
		instance = this;
		this.system = this.GetComponent<ParticleSystem>();
		this.RESOLUTION = new int[3] {(int) (GRID_RESOLUTION.x+2), (int) (GRID_RESOLUTION.y+2), (int) (GRID_RESOLUTION.z+2)};
		this.cells = new Cell[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		this.particles = new ParticleSystem.Particle[(RESOLUTION[0]) * (RESOLUTION[1]) * (RESOLUTION[2])];
		this.ORIGIN = this.transform.position;
		
		this.PARTICLE_SIZE = GRID_SIZE / Mathf.Max(GRID_RESOLUTION.x, GRID_RESOLUTION.y, GRID_RESOLUTION.z);
		InitializeGrid();

		if (UseTexture) {
			DIFFUSION_RATE /= 100;
			for (int x = 0; x < RESOLUTION[0]; x++) {
				for (int y = 0; y < RESOLUTION[1]; y++) {
					for (int z = 0; z < RESOLUTION[2]; z++) {
						Color pixelColor = ImportedTexture.GetPixel(x, z);
						cells[x, y, z].density = new float[3] { pixelColor.r, pixelColor.g, pixelColor.b };
						//cells[x, y, z].UpdateColor();
					}
				}
			}
		}

		if (ExternalForce) {
			ApplyExternalForce();
		}
	}

	void Update() {
		this.dt = Time.deltaTime;

		DetectInput();

		Diffuse();

		Advect();

		Project();

		Emit();

		UpdateParticleSystem();
	}


	public void ApplyExternalForce() {
		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					cells[i, j, k].velocity += externalForce;
				}
			}
		}
	}

	public void InitializeGrid() {
		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					Vector3 pos = ORIGIN + new Vector3(i, j, k) * PARTICLE_SIZE;
					cells[i, j, k] = new Cell();
					cells[i, j, k].position = pos;
					cells[i, j, k].index = new int[] {i, j, k};

					bool left = (i == 0);
					bool bottom = (j == 0);
					bool front = (k == 0);

					bool right = (i == RESOLUTION[0]-1);
					bool top = (j == RESOLUTION[1]-1);
					bool back = (k == RESOLUTION[2]-1);

					if (left || bottom || front || right || top || back) {
						cells[i, j, k].SetAsBoundary();
					}

					ParticleSystem.Particle p = new ParticleSystem.Particle();
					p.size = PARTICLE_SIZE;
					p.position = pos;
					cells[i, j, k].particle = p;
				}
			}
		}
	}

	public void UpdateParticleSystem() {
		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					UpdateParticle(i, j, k);
				}
			}
		}

		system.SetParticles(particles, particles.Length);
	}

	public void UpdateParticle(int i, int j, int k) {
		particles[(i) + (j) * (RESOLUTION[0]) + (k) * (RESOLUTION[1]) * (RESOLUTION[2])] = cells[i, j, k].particle;
		if (cells[i, j, k].isBoundary) {
			Debug.DrawLine(cells[i, j, k].position, cells[i, j, k].position + cells[i, j, k].velocity, Color.blue);
		} else {
			Debug.DrawLine(cells[i, j, k].position, cells[i, j, k].position + cells[i, j, k].velocity, Color.red);
		}
	}

	public ParticleSystem.Particle GetParticle(int i, int j, int k) {
		return cells[i, j, k].particle;
	}

	public void SetParticle(int i, int j, int k, ParticleSystem.Particle p) {
		cells[i, j, k].particle = p;
	}

	public void DetectInput() {
		if (Input.GetKeyUp("r")) {
			ChosenColor = 0;
		} else if (Input.GetKeyUp("g")) {
			ChosenColor = 1;
		} else if (Input.GetKeyUp("b")) {
			ChosenColor = 2;
		}

		for (int x = 0; x < RESOLUTION[0]; x++) {
			for (int y = 0; y < RESOLUTION[1]; y++) {
				for (int z = 0; z < RESOLUTION[2]; z++) {
					cells[x, y, z].prev_velocity = Vector3.zero;
					cells[x, y, z].prev_density = new float[3] { 0, 0, 0 };
					UpdateParticle(x, y, z);
				}
			}
		}

		int i = Mathf.RoundToInt((UserObject.transform.position.x - ORIGIN.x) / (PARTICLE_SIZE));
		int j = Mathf.RoundToInt((UserObject.transform.position.y - ORIGIN.y) / (PARTICLE_SIZE));
		int k = Mathf.RoundToInt((UserObject.transform.position.z - ORIGIN.z) / (PARTICLE_SIZE));
		//Debug.Log("" + i + " " + j + " " + k);
		if (i > 0 && i < RESOLUTION[0]-1 && j > 0 && j < RESOLUTION[1]-1 && k > 0 && k < RESOLUTION[2]-1) {
			Vector3 userVelocity = UserObject.GetComponent<CharacterController>().velocity;
			sources.Add(new Source(new int[] {i, j, k}, userVelocity, SourceStrength * userVelocity.magnitude, ChosenColor));
			//Debug.Log("within cell " + i + " " + j + " " + k);
			//Debug.Log(userVelocity);
		}
	}

	public void Emit() {
		/*foreach (Source s in sources) {
			int i = s.index[0];
			int j = s.index[1];
			int k = s.index[2];
			s.UpdateEmission	
		}*/

		List<Source> keep = new List<Source>();

		for (int a = 0; a < sources.Count; a++) {
			if (sources[a].isAlive()) {
				sources[a].UpdateEmission();
				keep.Add(sources[a]);
			}
		}

		sources = keep;
		//	Debug.Log(sources.Count);
	}

	public void Diffuse() {
		float a = dt * DIFFUSION_RATE * RESOLUTION[0] * RESOLUTION[1] * RESOLUTION[2];

		float[,,] d0 = new float[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		float[,,] d1 = new float[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		float[,,] d2 = new float[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		Vector3[,,] v0 = new Vector3[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];

		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					d0[i, j, k] = cells[i, j, k].density[0];
					d1[i, j, k] = cells[i, j, k].density[1];
					d2[i, j, k] = cells[i, j, k].density[2];
					v0[i, j, k] = cells[i, j, k].velocity;
				}
			}
		}

		// Gauss-Seidel relaxation
		for (int r = 0; r < relaxationIterations; r++) {
			for (int i = 1; i < RESOLUTION[0]-1; i++) {
				for (int j = 1; j < RESOLUTION[1]-1; j++) {
					for (int k = 1; k < RESOLUTION[2]-1; k++) {
						float newDensity0 = (d0[i, j, k] + a * (cells[i+1, j, k].density[0] + cells[i-1, j, k].density[0] + cells[i, j+1, k].density[0] + cells[i, j-1, k].density[0] + cells[i, j, k+1].density[0] + cells[i, j, k-1].density[0])) / (1+6*a);
						float newDensity1 = (d1[i, j, k] + a * (cells[i+1, j, k].density[1] + cells[i-1, j, k].density[1] + cells[i, j+1, k].density[1] + cells[i, j-1, k].density[1] + cells[i, j, k+1].density[1] + cells[i, j, k-1].density[1])) / (1+6*a);
						float newDensity2 = (d2[i, j, k] + a * (cells[i+1, j, k].density[2] + cells[i-1, j, k].density[2] + cells[i, j+1, k].density[2] + cells[i, j-1, k].density[2] + cells[i, j, k+1].density[2] + cells[i, j, k-1].density[2])) / (1+6*a);
						Vector3 newVelocity = (v0[i, j, k] + a * (cells[i+1, j, k].velocity + cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity + cells[i, j, k-1].velocity)) / (1+6*a);

						//f (r == relaxationIterations-1) {
							cells[i, j, k].SetDensity(newDensity0, 0);
							cells[i, j, k].SetDensity(newDensity1, 1);
							cells[i, j, k].SetDensity(newDensity2, 2);
							cells[i, j, k].SetVelocity(newVelocity);
						//} else {
						//	cells[i, j, k].density = newDensity;
						//	cells[i, j, k].velocity = newVelocity;
						//}				
					}
				}
			}
			SetDensityBoundaries(0);
			SetDensityBoundaries(1);
			SetDensityBoundaries(2);
			SetVelocityBoundaries();
		}
	}

	public void Advect() {
		float dt0i = dt * RESOLUTION[0];
		float dt0j = dt * RESOLUTION[1];
		float dt0k = dt * RESOLUTION[2];

		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					float x = Mathf.Clamp(i - dt0i * cells[i, j, k].velocity.x, 0.5f, RESOLUTION[0]-1.5f);
					float y = Mathf.Clamp(j - dt0j * cells[i, j, k].velocity.y, 0.5f, RESOLUTION[1]-1.5f);
					float z = Mathf.Clamp(k - dt0k * cells[i, j, k].velocity.z, 0.5f, RESOLUTION[2]-1.5f);

					int i0 = (int) (x); 	int i1 = i0+1;
					int j0 = (int) (y); 	int j1 = j0+1;
					int k0 = (int) (z); 	int k1 = k0+1;

					float s1 = x - i0;	float s0 = 1-s1;
					float t1 = y - j0;	float t0 = 1-t1;
					float u1 = z - k0;	float u0 = 1-u1;

					//Debug.Log(i0 + " " + j0 + " " + k0);
					float r000 = cells[i0, j0, k0].prev_density[0];
					float r100 = cells[i1, j0, k0].prev_density[0];
					float r010 = cells[i0, j1, k0].prev_density[0];
					float r110 = cells[i1, j1, k0].prev_density[0];
					float r001 = cells[i0, j0, k1].prev_density[0];
					float r101 = cells[i1, j0, k1].prev_density[0];
					float r011 = cells[i0, j1, k1].prev_density[0];
					float r111 = cells[i1, j1, k1].prev_density[0];

					float g000 = cells[i0, j0, k0].prev_density[1];
					float g100 = cells[i1, j0, k0].prev_density[1];
					float g010 = cells[i0, j1, k0].prev_density[1];
					float g110 = cells[i1, j1, k0].prev_density[1];
					float g001 = cells[i0, j0, k1].prev_density[1];
					float g101 = cells[i1, j0, k1].prev_density[1];
					float g011 = cells[i0, j1, k1].prev_density[1];
					float g111 = cells[i1, j1, k1].prev_density[1];

					float b000 = cells[i0, j0, k0].prev_density[2];
					float b100 = cells[i1, j0, k0].prev_density[2];
					float b010 = cells[i0, j1, k0].prev_density[2];
					float b110 = cells[i1, j1, k0].prev_density[2];
					float b001 = cells[i0, j0, k1].prev_density[2];
					float b101 = cells[i1, j0, k1].prev_density[2];
					float b011 = cells[i0, j1, k1].prev_density[2];
					float b111 = cells[i1, j1, k1].prev_density[2];

					Vector3 v000 = cells[i0, j0, k0].prev_velocity;
					Vector3 v100 = cells[i1, j0, k0].prev_velocity;
					Vector3 v010 = cells[i0, j1, k0].prev_velocity;
					Vector3 v110 = cells[i1, j1, k0].prev_velocity;
					Vector3 v001 = cells[i0, j0, k1].prev_velocity;
					Vector3 v101 = cells[i1, j0, k1].prev_velocity;
					Vector3 v011 = cells[i0, j1, k1].prev_velocity;
					Vector3 v111 = cells[i1, j1, k1].prev_velocity;

					float newDensity0 = s0 * (t0 * (u0 * r000 + u1 * r001) + t1 * (u0 * r010 + u1 * r011)) + s1 * (t0 * (u0 * r100 + u1 * r101) + t1 * (u0 * r110 + u1 * r111));
				 	cells[i, j, k].SetDensity(newDensity0, 0);

				 	float newDensity1 = s0 * (t0 * (u0 * g000 + u1 * g001) + t1 * (u0 * g010 + u1 * g011)) + s1 * (t0 * (u0 * g100 + u1 * g101) + t1 * (u0 * g110 + u1 * g111));
				 	cells[i, j, k].SetDensity(newDensity1, 1);

				 	float newDensity2 = s0 * (t0 * (u0 * b000 + u1 * b001) + t1 * (u0 * b010 + u1 * b011)) + s1 * (t0 * (u0 * b100 + u1 * b101) + t1 * (u0 * b110 + u1 * b111));
				 	cells[i, j, k].SetDensity(newDensity2, 2);

				 	Vector3 newVelocity = s0 * (t0 * (u0 * v000 + u1 * v001) + t1 * (u0 * v010 + u1 * v011)) + s1 * (t0 * (u0 * v100 + u1 * v101) + t1 * (u0 * v110 + u1 * v111));
					cells[i, j, k].SetVelocity(newVelocity);
				}
			}
		}

		SetDensityBoundaries(0);
		SetDensityBoundaries(1);
		SetDensityBoundaries(2);
	}

	public void Project() {
		Vector3 h = new Vector3(1f/RESOLUTION[0], 1f/RESOLUTION[1], 1f/RESOLUTION[2]);
		Vector3[,,] f0 = new Vector3[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];

		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					f0[i, j, k] = -0.5f * (h.x * (cells[i+1, j, k].velocity - cells[i-1, j, k].velocity) + h.y * (cells[i, j+1, k].velocity - cells[i, j-1, k].velocity) + h.z * (cells[i, j, k+1].velocity - cells[i, j, k-1].velocity));
				}
			}
		}

		// Gauss-Seidel relaxation
		for (int r = 0; r < relaxationIterations; r++) {
			for (int i = 1; i < RESOLUTION[0]-1; i++) {
				for (int j = 1; j < RESOLUTION[1]-1; j++) {
					for (int k = 1; k < RESOLUTION[2]-1; k++) {
						cells[i, j, k].pressure = (f0[i, j, k] + cells[i-1, j, k].pressure + cells[i+1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k-1].pressure + cells[i, j, k+1].pressure) / 6f;
					}
				}
			}
			SetPressureBoundaries();
		}



		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					Vector3 newVelocity = cells[i, j, k].velocity;
					newVelocity.x -= 0.5f * (cells[i+1, j, k].pressure.x - cells[i-1, j, k].pressure.x) / h.x;
					newVelocity.y -= 0.5f * (cells[i, j+1, k].pressure.y - cells[i, j-1, k].pressure.y) / h.y;
					newVelocity.z -= 0.5f * (cells[i, j, k+1].pressure.z - cells[i, j, k-1].pressure.z) / h.z;

					cells[i, j, k].SetVelocity(newVelocity);
				}
			}
		}
	}

	/*public void SetBoundaries() {
		for (int i = 0; i < RESOLUTION[0]; i += RESOLUTION[0]-1) {
			for (int j = 0; j < RESOLUTION[1]; j += RESOLUTION[1]-1) {
				for (int k = 0; k < RESOLUTION[2]; k += RESOLUTION[2]-1) {
					bool left = (i == 0);
					bool bottom = (j == 0);
					bool front = (k == 0);

					bool right = (i == RESOLUTION[0]-1);
					bool top = (j == RESOLUTION[1]-1);
					bool back = (k == RESOLUTION[2]-1);
					if (left || right || bottom || front || top || back) {
						Vector3 newVelocity, newPressure;
						float newDensity;

						if (left && bottom && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j+1, k].density + cells[i, j, k+1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i+1, j, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k+1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && bottom && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j+1, k].density + cells[i, j, k-1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i+1, j, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k-1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && top && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j-1, k].density + cells[i, j, k+1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i+1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j, k+1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && top && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j-1, k].density + cells[i, j, k-1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i+1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j, k-1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && bottom && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j+1, k].density + cells[i, j, k+1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i-1, j, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k+1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && bottom && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j+1, k].density + cells[i, j, k-1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i-1, j, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k-1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && top && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j-1, k].density + cells[i, j, k+1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i-1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j, k+1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && top && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j-1, k].density + cells[i, j, k-1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i-1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j, k-1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && bottom) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j+1, k].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i+1, j, k].pressure + cells[i, j+1, k].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && top) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j-1, k].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i+1, j, k].pressure + cells[i, j-1, k].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j, k+1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i+1, j, k].pressure + cells[i, j, k+1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j, k-1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i+1, j, k].pressure + cells[i, j, k-1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && bottom) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j+1, k].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i-1, j, k].pressure + cells[i, j+1, k].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && top) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j-1, k].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i-1, j, k].pressure + cells[i, j-1, k].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j, k+1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i-1, j, k].pressure + cells[i, j, k+1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j, k-1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i-1, j, k].pressure + cells[i, j, k-1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (bottom && front) {
							newVelocity = (cells[i, j+1, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i, j+1, k].density + cells[i, j, k+1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i, j+1, k].pressure + cells[i, j, k+1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (bottom && back) {
							newVelocity = (cells[i, j+1, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i, j+1, k].density + cells[i, j, k-1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i, j+1, k].pressure + cells[i, j, k-1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (top && front) {
							newVelocity = (cells[i, j-1, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i, j-1, k].density + cells[i, j, k+1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i, j-1, k].pressure + cells[i, j, k+1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (top && back) {
							newVelocity = (cells[i, j-1, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i, j-1, k].density + cells[i, j, k-1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
							newPressure = (cells[i, j-1, k].pressure + cells[i, j, k-1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (left) {
							newVelocity.x = -1 *cells[i+1, j, k].velocity.x;
							newVelocity.y = 	cells[i+1, j, k].velocity.y;
							newVelocity.z = 	cells[i+1, j, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i+1, j, k].density);
							cells[i, j, k].pressure = cells[i+1, j, k].pressure;
						} else if (right) {
							newVelocity.x = -1 *cells[i-1, j, k].velocity.x;
							newVelocity.y = 	cells[i-1, j, k].velocity.y;
							newVelocity.z = 	cells[i-1, j, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i-1, j, k].density);
							cells[i, j, k].pressure = cells[i-1, j, k].pressure;
						} else if (bottom) {
							newVelocity.x = 	cells[i, j+1, k].velocity.x;
							newVelocity.y = -1 *cells[i, j+1, k].velocity.y;
							newVelocity.z = 	cells[i, j+1, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i, j+1, k].density);
							cells[i, j, k].pressure = cells[i, j+1, k].pressure;
						} else if (top) {
							newVelocity.x = 	cells[i, j-1, k].velocity.x;
							newVelocity.y = -1 *cells[i, j-1, k].velocity.y;
							newVelocity.z = 	cells[i, j-1, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i, j-1, k].density);
							cells[i, j, k].pressure = cells[i, j-1, k].pressure;
						} else if (front) {
							newVelocity.x = 	cells[i, j, k+1].velocity.x;
							newVelocity.y = 	cells[i, j, k+1].velocity.y;
							newVelocity.z = -1 *cells[i, j, k+1].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i, j, k+1].density);
							cells[i, j, k].pressure = cells[i, j, k+1].pressure;
						} else if (back) {
							newVelocity.x = 	cells[i, j, k-1].velocity.x;
							newVelocity.y =		cells[i, j, k-1].velocity.y;
							newVelocity.z = -1 *cells[i, j, k-1].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i, j, k-1].density);
							cells[i, j, k].pressure = cells[i, j, k-1].pressure;
						}
						Debug.DrawLine(cells[i, j, k].position, cells[i, j, k].position + cells[i, j, k].velocity, Color.blue);
					}
				}
			}
		}
	}*/

	public void SetDensityBoundaries(int ci) {
		for (int i = 0; i < RESOLUTION[0]; i += RESOLUTION[0]-1) {
			for (int j = 0; j < RESOLUTION[1]; j += RESOLUTION[1]-1) {
				for (int k = 0; k < RESOLUTION[2]; k += RESOLUTION[2]-1) {
					bool left = (i == 0);
					bool bottom = (j == 0);
					bool front = (k == 0);

					bool right = (i == RESOLUTION[0]-1);
					bool top = (j == RESOLUTION[1]-1);
					bool back = (k == RESOLUTION[2]-1);
					if (left || right || bottom || front || top || back) {
						float newDensity;

						if (left && bottom && front) {
							newDensity = (cells[i+1, j, k].density[ci] + cells[i, j+1, k].density[ci] + cells[i, j, k+1].density[ci])/3f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (left && bottom && back) {
							newDensity = (cells[i+1, j, k].density[ci] + cells[i, j+1, k].density[ci] + cells[i, j, k-1].density[ci])/3f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (left && top && front) {
							newDensity = (cells[i+1, j, k].density[ci] + cells[i, j-1, k].density[ci] + cells[i, j, k+1].density[ci])/3f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (left && top && back) {
							newDensity = (cells[i+1, j, k].density[ci] + cells[i, j-1, k].density[ci] + cells[i, j, k-1].density[ci])/3f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (right && bottom && front) {
							newDensity = (cells[i-1, j, k].density[ci] + cells[i, j+1, k].density[ci] + cells[i, j, k+1].density[ci])/3f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (right && bottom && back) {
							newDensity = (cells[i-1, j, k].density[ci] + cells[i, j+1, k].density[ci] + cells[i, j, k-1].density[ci])/3f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (right && top && front) {
							newDensity = (cells[i-1, j, k].density[ci] + cells[i, j-1, k].density[ci] + cells[i, j, k+1].density[ci])/3f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (right && top && back) {
							newDensity = (cells[i-1, j, k].density[ci] + cells[i, j-1, k].density[ci] + cells[i, j, k-1].density[ci])/3f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (left && bottom) {
							newDensity = (cells[i+1, j, k].density[ci] + cells[i, j+1, k].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (left && top) {
							newDensity = (cells[i+1, j, k].density[ci] + cells[i, j-1, k].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (left && front) {
							newDensity = (cells[i+1, j, k].density[ci] + cells[i, j, k+1].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (left && back) {
							newDensity = (cells[i+1, j, k].density[ci] + cells[i, j, k-1].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (right && bottom) {
							newDensity = (cells[i-1, j, k].density[ci] + cells[i, j+1, k].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (right && top) {
							newDensity = (cells[i-1, j, k].density[ci] + cells[i, j-1, k].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (right && front) {
							newDensity = (cells[i-1, j, k].density[ci] + cells[i, j, k+1].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (right && back) {
							newDensity = (cells[i-1, j, k].density[ci] + cells[i, j, k-1].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (bottom && front) {
							newDensity = (cells[i, j+1, k].density[ci] + cells[i, j, k+1].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (bottom && back) {
							newDensity = (cells[i, j+1, k].density[ci] + cells[i, j, k-1].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (top && front) {
							newDensity = (cells[i, j-1, k].density[ci] + cells[i, j, k+1].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (top && back) {
							newDensity = (cells[i, j-1, k].density[ci] + cells[i, j, k-1].density[ci])/2f;
							cells[i, j, k].SetDensity(newDensity, ci);
						} else if (left) {
							cells[i, j, k].SetDensity(cells[i+1, j, k].density[ci], ci);
						} else if (right) {
							cells[i, j, k].SetDensity(cells[i-1, j, k].density[ci], ci);
						} else if (bottom) {
							cells[i, j, k].SetDensity(cells[i, j+1, k].density[ci], ci);
						} else if (top) {
							cells[i, j, k].SetDensity(cells[i, j-1, k].density[ci], ci);
						} else if (front) {
							cells[i, j, k].SetDensity(cells[i, j, k+1].density[ci], ci);
						} else if (back) {
							cells[i, j, k].SetDensity(cells[i, j, k-1].density[ci], ci);
						}
					}
				}
			}
		}
	}

	public void SetVelocityBoundaries() {
		for (int i = 0; i < RESOLUTION[0]; i += RESOLUTION[0]-1) {
			for (int j = 0; j < RESOLUTION[1]; j += RESOLUTION[1]-1) {
				for (int k = 0; k < RESOLUTION[2]; k += RESOLUTION[2]-1) {
					bool left = (i == 0);
					bool bottom = (j == 0);
					bool front = (k == 0);

					bool right = (i == RESOLUTION[0]-1);
					bool top = (j == RESOLUTION[1]-1);
					bool back = (k == RESOLUTION[2]-1);
					if (left || right || bottom || front || top || back) {
						Vector3 newVelocity;

						if (left && bottom && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (left && bottom && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (left && top && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (left && top && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right && bottom && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right && bottom && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right && top && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right && top && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (left && bottom) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (left && top) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (left && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (left && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right && bottom) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right && top) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (bottom && front) {
							newVelocity = (cells[i, j+1, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (bottom && back) {
							newVelocity = (cells[i, j+1, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (top && front) {
							newVelocity = (cells[i, j-1, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (top && back) {
							newVelocity = (cells[i, j-1, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (left) {
							newVelocity.x = -1 *cells[i+1, j, k].velocity.x;
							newVelocity.y = 	cells[i+1, j, k].velocity.y;
							newVelocity.z = 	cells[i+1, j, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (right) {
							newVelocity.x = -1 *cells[i-1, j, k].velocity.x;
							newVelocity.y = 	cells[i-1, j, k].velocity.y;
							newVelocity.z = 	cells[i-1, j, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (bottom) {
							newVelocity.x = 	cells[i, j+1, k].velocity.x;
							newVelocity.y = -1 *cells[i, j+1, k].velocity.y;
							newVelocity.z = 	cells[i, j+1, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (top) {
							newVelocity.x = 	cells[i, j-1, k].velocity.x;
							newVelocity.y = -1 *cells[i, j-1, k].velocity.y;
							newVelocity.z = 	cells[i, j-1, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (front) {
							newVelocity.x = 	cells[i, j, k+1].velocity.x;
							newVelocity.y = 	cells[i, j, k+1].velocity.y;
							newVelocity.z = -1 *cells[i, j, k+1].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
						} else if (back) {
							newVelocity.x = 	cells[i, j, k-1].velocity.x;
							newVelocity.y =		cells[i, j, k-1].velocity.y;
							newVelocity.z = -1 *cells[i, j, k-1].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
						}
						Debug.DrawLine(cells[i, j, k].position, cells[i, j, k].position + cells[i, j, k].velocity, Color.blue);
					}
				}
			}
		}
	}

	public void SetPressureBoundaries() {
		for (int i = 0; i < RESOLUTION[0]; i += RESOLUTION[0]-1) {
			for (int j = 0; j < RESOLUTION[1]; j += RESOLUTION[1]-1) {
				for (int k = 0; k < RESOLUTION[2]; k += RESOLUTION[2]-1) {
					bool left = (i == 0);
					bool bottom = (j == 0);
					bool front = (k == 0);

					bool right = (i == RESOLUTION[0]-1);
					bool top = (j == RESOLUTION[1]-1);
					bool back = (k == RESOLUTION[2]-1);
					if (left || right || bottom || front || top || back) {
						Vector3 newPressure;

						if (left && bottom && front) {
							newPressure = (cells[i+1, j, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k+1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && bottom && back) {
							newPressure = (cells[i+1, j, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k-1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && top && front) {
							newPressure = (cells[i+1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j, k+1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && top && back) {
							newPressure = (cells[i+1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j, k-1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && bottom && front) {
							newPressure = (cells[i-1, j, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k+1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && bottom && back) {
							newPressure = (cells[i-1, j, k].pressure + cells[i, j+1, k].pressure + cells[i, j, k-1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && top && front) {
							newPressure = (cells[i-1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j, k+1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && top && back) {
							newPressure = (cells[i-1, j, k].pressure + cells[i, j-1, k].pressure + cells[i, j, k-1].pressure)/3f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && bottom) {
							newPressure = (cells[i+1, j, k].pressure + cells[i, j+1, k].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && top) {
							newPressure = (cells[i+1, j, k].pressure + cells[i, j-1, k].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && front) {
							newPressure = (cells[i+1, j, k].pressure + cells[i, j, k+1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (left && back) {
							newPressure = (cells[i+1, j, k].pressure + cells[i, j, k-1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && bottom) {
							newPressure = (cells[i-1, j, k].pressure + cells[i, j+1, k].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && top) {
							newPressure = (cells[i-1, j, k].pressure + cells[i, j-1, k].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && front) {
							newPressure = (cells[i-1, j, k].pressure + cells[i, j, k+1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (right && back) {
							newPressure = (cells[i-1, j, k].pressure + cells[i, j, k-1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (bottom && front) {
							newPressure = (cells[i, j+1, k].pressure + cells[i, j, k+1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (bottom && back) {
							newPressure = (cells[i, j+1, k].pressure + cells[i, j, k-1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (top && front) {
							newPressure = (cells[i, j-1, k].pressure + cells[i, j, k+1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (top && back) {
							newPressure = (cells[i, j-1, k].pressure + cells[i, j, k-1].pressure)/2f;
							cells[i, j, k].pressure = newPressure;
						} else if (left) {
							cells[i, j, k].pressure = cells[i+1, j, k].pressure;
						} else if (right) {
							cells[i, j, k].pressure = cells[i-1, j, k].pressure;
						} else if (bottom) {
							cells[i, j, k].pressure = cells[i, j+1, k].pressure;
						} else if (top) {
							cells[i, j, k].pressure = cells[i, j-1, k].pressure;
						} else if (front) {
							cells[i, j, k].pressure = cells[i, j, k+1].pressure;
						} else if (back) {
							cells[i, j, k].pressure = cells[i, j, k-1].pressure;
						}
					}
				}
			}
		}
	}
}