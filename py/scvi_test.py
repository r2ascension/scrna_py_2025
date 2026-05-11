#!/usr/bin/env python3
"""
scvi_integration_gpu_optimized.py

全面优化的scVI批次整合脚本
- GPU加速训练
- 高变基因筛选
- 批次整合效果评估
- 多分辨率聚类
- 增强可视化

作者:临床-生信团队
日期:2025-10-30
优化版本:v2.0
"""

import sys
import os
from pathlib import Path
import warnings
import numpy as np
warnings.filterwarnings('ignore')

# ==================== 配置基础 ====================

# ========== 输入输出配置 ==========
INPUT_H5AD_PATH = "/home/h2048/data/R/1029/seurat_without_partial_dual_contamination_doublets_2_1029.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1029/scvi_output_optimized"
OVERWRITE_EXISTING = True

# ========== GPU配置 ==========
USE_GPU = True  # 是否使用GPU(如果可用)
GPU_DEVICE = 0  # 使用哪个GPU(0, 1, 2...)或"auto"自动选择

# ========== 模型加载配置 ==========
# 是否加载已有的scVI模型(而不是重新训练)
LOAD_EXISTING_MODEL = False  # 设为True可跳过训练,直接加载模型
EXISTING_MODEL_PATH = "/home/h2048/data/py/1029/scvi_output_optimized/scvi_model"  # 已训练模型路径

# 是否直接读取已整合的h5ad数据(跳过所有整合步骤,直接进行下游分析)
LOAD_INTEGRATED_DATA = False  # 设为True可跳过整合,直接分析
INTEGRATED_H5AD_PATH = "/home/h2048/data/py/1029/scvi_output_optimized/adata_scvi_integrated.h5ad"  # 已整合数据路径

# ========== scVI整合配置 ==========
BATCH_KEY = "dataset"  # 批次变量名称

# scVI模型参数
SCVI_N_LAYERS = 4  # 编码器/解码器层数（1-5，推荐2-4）
SCVI_N_LATENT = 75  # latent空间维度（10-100，推荐30-50）
SCVI_GENE_LIKELIHOOD = "nb"  # 基因似然模型（"nb"或"zinb"）
SCVI_DROPOUT_RATE = 0.1  # dropout率（0.0-0.5，推荐0.1）
SCVI_DISPERSION = "gene"  # dispersion参数（"gene"或"gene-batch"）

# scVI模型参数对比（可选，用于找最佳配置）
COMPARE_SCVI_PARAMS = False  # 是否训练多组scVI参数进行对比
SCVI_PARAM_GRID = [
    {"n_layers": 2, "n_latent": 30, "name": "light"},      # 轻量模型
    {"n_layers": 3, "n_latent": 50, "name": "balanced"},   # 平衡模型（当前）
    {"n_layers": 3, "n_latent": 75, "name": "rich"},       # 丰富latent
    {"n_layers": 4, "n_latent": 50, "name": "deep"},       # 深层模型
]
SCVI_MAX_EPOCHS = 500  # 明确指定最大轮数

# 训练参数(已优化)
TRAIN_BATCH_SIZE = 256  # 增大batch size(GPU显存足够时)
TRAIN_LR = 1e-3  # 学习率
EARLY_STOPPING = True  # 使用early stopping
EARLY_STOPPING_PATIENCE = 50  # 增大耐心值

# ========== 高变基因筛选 ==========
USE_HVG = True  # 是否筛选高变基因
N_TOP_GENES = 3000  # 保留的高变基因数量
HVG_FLAVOR = "seurat_v3"  # 高变基因筛选方法

# ========== 连续协变量 ==========
CONTINUOUS_COVARIATES = [
    # "percent.mt",
    # "percent.rp",
]

# ========== 基因过滤配置 ==========
MIN_CELLS_PER_GENE = 3

# ========== 降维和可视化配置 ==========
RUN_UMAP = True
# UMAP参数优化（针对批次整合后的可视化）
UMAP_MIN_DIST = 0.5  # 降低以获得更紧凑的嵌入（原0.3）
UMAP_N_NEIGHBORS = 50  # 增加以捕获更多全局结构（原30）
UMAP_SPREAD = 1.0  # 控制嵌入点的分散程度
UMAP_N_COMPONENTS = 2  # 降维维度（2D可视化）
# neighbors计算使用的latent维度
NEIGHBORS_N_PCS = None  # None则使用全部latent维度

# 多参数UMAP比较（可选，用于参数调优）
COMPARE_UMAP_PARAMS = False  # 是否生成多组参数的UMAP比较图
UMAP_PARAM_GRID = [
    {"n_neighbors": 30, "min_dist": 0.1, "name": "tight_local"},
    {"n_neighbors": 50, "min_dist": 0.1, "name": "tight_global"},
    {"n_neighbors": 30, "min_dist": 0.3, "name": "loose_local"},
    {"n_neighbors": 50, "min_dist": 0.3, "name": "loose_global"},
]

# ========== 聚类配置(多分辨率)==========
RUN_CLUSTERING = True
LEIDEN_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8]  # 多个分辨率
DEFAULT_RESOLUTION = 0.4  # 默认使用的分辨率

# ========== 批次整合评估 ==========
EVALUATE_INTEGRATION = True  # 是否评估批次整合效果
ANALYZE_UNINTEGRATED = False  # 是否分析未整合数据(评估批次效应)

# ========== 可视化变量 ==========
VISUALIZATION_VARS = [
    "dataset",
    "Annotation",
    "tissue_sampling_method",
]

# ========== 增强可视化 ==========
GENERATE_FACET_PLOTS = True  # 是否生成分批次显示图

VERBOSE = True

# ==================== 主程序 ====================

def log_msg(msg):
    """打印日志"""
    if VERBOSE:
        print(msg)


def setup_gpu():
    """配置GPU环境"""
    import torch
    
    log_msg("\n" + "="*70)
    log_msg("🔧 配置GPU环境")
    log_msg("="*70)
    
    # 检查CUDA是否可用
    if not torch.cuda.is_available():
        log_msg("\n⚠️  CUDA不可用")
        log_msg(f"   PyTorch版本: {torch.__version__}")
        log_msg(f"   CUDA编译版本: {torch.version.cuda}")
        
        if USE_GPU:
            log_msg("\n❌ 配置要求使用GPU,但CUDA不可用")
            log_msg("\n💡 可能的原因:")
            log_msg("   1. 安装了CPU版本的PyTorch")
            log_msg("   2. CUDA驱动未正确安装")
            log_msg("   3. PyTorch CUDA版本与系统不匹配")
            log_msg("\n💡 解决方案:")
            log_msg("   运行诊断脚本: python gpu_diagnostics.py")
            log_msg("   或重新安装PyTorch CUDA版本")
            sys.exit(1)
        else:
            log_msg("   将使用CPU运行训练(速度较慢)")
            return None
    
    # GPU可用
    n_gpus = torch.cuda.device_count()
    log_msg(f"\n✅ 检测到 {n_gpus} 个GPU:")
    
    for i in range(n_gpus):
        gpu_name = torch.cuda.get_device_name(i)
        gpu_memory = torch.cuda.get_device_properties(i).total_memory / 1e9
        log_msg(f"   GPU {i}: {gpu_name} ({gpu_memory:.1f} GB)")
    
    # 选择GPU设备
    if GPU_DEVICE == "auto":
        # 自动选择显存最多的GPU
        gpu_id = 0
        max_memory = 0
        for i in range(n_gpus):
            memory = torch.cuda.get_device_properties(i).total_memory
            if memory > max_memory:
                max_memory = memory
                gpu_id = i
        log_msg(f"\n🎯 自动选择GPU {gpu_id}")
    elif isinstance(GPU_DEVICE, str) and "," in GPU_DEVICE:
        # 多GPU
        gpu_ids = [int(x.strip()) for x in GPU_DEVICE.split(",")]
        log_msg(f"\n🎯 使用多GPU: {gpu_ids}")
        gpu_id = gpu_ids[0]  # 主GPU
    else:
        gpu_id = int(GPU_DEVICE)
        log_msg(f"\n🎯 使用GPU {gpu_id}")
    
    # 设置默认GPU
    torch.cuda.set_device(gpu_id)
    
    # 显示当前GPU信息
    current_device = torch.cuda.current_device()
    log_msg(f"\n📊 当前GPU信息:")
    log_msg(f"   设备ID: {current_device}")
    log_msg(f"   设备名称: {torch.cuda.get_device_name(current_device)}")
    log_msg(f"   总显存: {torch.cuda.get_device_properties(current_device).total_memory / 1e9:.2f} GB")
    
    # 测试GPU
    try:
        test_tensor = torch.zeros(1).cuda()
        del test_tensor
        torch.cuda.empty_cache()
        log_msg(f"   ✅ GPU测试成功")
    except Exception as e:
        log_msg(f"   ❌ GPU测试失败: {e}")
        if USE_GPU:
            sys.exit(1)
    
    return gpu_id


def load_anndata(h5ad_path):
    """读取h5ad数据"""
    import scanpy as sc
    
    log_msg("\n" + "="*70)
    log_msg("📂 Step 1: Loading data")
    log_msg("="*70)
    log_msg(f"\n📁 Reading: {h5ad_path}")
    
    if not Path(h5ad_path).exists():
        raise FileNotFoundError(f"File not found: {h5ad_path}")
    
    # 读取数据
    adata = sc.read_h5ad(h5ad_path)
    
    log_msg(f"✅ Data loaded successfully")
    log_msg(f"   Cells: {adata.n_obs:,}")
    log_msg(f"   Genes: {adata.n_vars:,}")
    
    # 显示obs列
    log_msg(f"\n📋 Available metadata columns:")
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        log_msg(f"   - {col}: {n_unique} unique values")
    
    # 检查batch key
    if BATCH_KEY not in adata.obs.columns:
        raise ValueError(f"Batch key '{BATCH_KEY}' not found in adata.obs")
    
    # 显示batch分布
    log_msg(f"\n📊 Batch distribution (key: {BATCH_KEY}):")
    batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        log_msg(f"   {batch}: {count:,} cells ({pct:.1f}%)")
    
    return adata


def preprocess_for_scvi(adata):
    """预处理数据用于scVI"""
    import scanpy as sc
    
    log_msg("\n" + "="*70)
    log_msg("🔬 Step 2: Preprocessing for scVI")
    log_msg("="*70)
    
    # 复制数据避免修改原始数据
    adata = adata.copy()
    
    # 基因过滤
    log_msg(f"\n🧬 Gene filtering (min_cells={MIN_CELLS_PER_GENE})...")
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    n_genes_after = adata.n_vars
    log_msg(f"   Genes before: {n_genes_before:,}")
    log_msg(f"   Genes after: {n_genes_after:,}")
    log_msg(f"   Removed: {n_genes_before - n_genes_after:,}")
    
    # 高变基因筛选
    if USE_HVG:
        log_msg(f"\n🎯 Selecting highly variable genes (n={N_TOP_GENES})...")
        log_msg(f"   Method: {HVG_FLAVOR}")
        
        # 保存raw counts
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()
        
        # 计算高变基因
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=N_TOP_GENES,
            flavor=HVG_FLAVOR,
            batch_key=BATCH_KEY,
            subset=False  # 先不subset,保留所有基因信息
        )
        
        n_hvg = adata.var['highly_variable'].sum()
        log_msg(f"   ✅ Selected {n_hvg} highly variable genes")
        
        # 子集到高变基因
        adata = adata[:, adata.var['highly_variable']].copy()
        log_msg(f"   Final genes: {adata.n_vars:,}")
    
    # 确保使用raw counts
    if 'counts' in adata.layers:
        log_msg("\n📊 Using raw counts from layers['counts']")
        adata.X = adata.layers['counts'].copy()
    else:
        log_msg("\n⚠️  No 'counts' layer found, using X as counts")
    
    log_msg("\n✅ Preprocessing completed")
    log_msg(f"   Final dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    return adata


def run_scvi_integration_gpu(adata, output_dir, gpu_id=None):
    """运行scVI批次整合(GPU加速版)"""
    import scanpy as sc
    import scvi
    import torch
    
    log_msg("\n" + "="*70)
    log_msg("🚀 Step 3: scVI Integration (GPU Accelerated)")
    log_msg("="*70)
    
    # 预处理
    adata = preprocess_for_scvi(adata)
    
    # 配置scVI
    log_msg("\n⚙️  Configuring scVI model...")
    log_msg(f"   n_layers: {SCVI_N_LAYERS}")
    log_msg(f"   n_latent: {SCVI_N_LATENT}")
    log_msg(f"   gene_likelihood: {SCVI_GENE_LIKELIHOOD}")
    log_msg(f"   dropout_rate: {SCVI_DROPOUT_RATE}")
    log_msg(f"   dispersion: {SCVI_DISPERSION}")
    log_msg(f"   batch_key: {BATCH_KEY}")
    log_msg(f"   continuous_covariates: {CONTINUOUS_COVARIATES if CONTINUOUS_COVARIATES else 'None'}")
    
    # 设置scVI的AnnData
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=BATCH_KEY,
        continuous_covariate_keys=CONTINUOUS_COVARIATES if CONTINUOUS_COVARIATES else None
    )
    
    # 模型路径
    model_save_path = output_dir / "scvi_model"
    
    # 决定是加载现有模型还是训练新模型
    # 使用局部变量避免UnboundLocalError
    should_load_model = LOAD_EXISTING_MODEL and model_save_path.exists()
    model_loaded = False
    
    if should_load_model:
        log_msg(f"\n📥 Loading existing scVI model from: {model_save_path}")
        try:
            model = scvi.model.SCVI.load(model_save_path, adata=adata)
            log_msg("   ✅ Model loaded successfully")
            model_loaded = True
        except Exception as e:
            log_msg(f"   ❌ Failed to load model: {e}")
            log_msg("   Will train a new model instead...")
            model_loaded = False
    
    if not model_loaded:
        # 创建模型
        log_msg("\n🏗️  Creating scVI model...")
        model = scvi.model.SCVI(
            adata,
            n_layers=SCVI_N_LAYERS,
            n_latent=SCVI_N_LATENT,
            gene_likelihood=SCVI_GENE_LIKELIHOOD,
            dropout_rate=SCVI_DROPOUT_RATE,
            dispersion=SCVI_DISPERSION,
        )
        
        # GPU设置
        if USE_GPU and gpu_id is not None:
            log_msg(f"\n🎮 Using GPU {gpu_id}: {torch.cuda.get_device_name(gpu_id)}")
            accelerator = "gpu"
            devices = [gpu_id]
        else:
            log_msg("\n💻 Using CPU (this will be slower)")
            accelerator = "cpu"
            devices = "auto"
        
        # 训练参数
        log_msg("\n🏋️  Training parameters:")
        log_msg(f"   max_epochs: {SCVI_MAX_EPOCHS}")
        log_msg(f"   batch_size: {TRAIN_BATCH_SIZE}")
        log_msg(f"   learning_rate: {TRAIN_LR}")
        log_msg(f"   early_stopping: {EARLY_STOPPING}")
        if EARLY_STOPPING:
            log_msg(f"   patience: {EARLY_STOPPING_PATIENCE}")
        
        # 训练模型
        log_msg("\n🚂 Training scVI model...")
        log_msg("   This may take several minutes depending on data size...")
        
        try:
            model.train(
                max_epochs=SCVI_MAX_EPOCHS,
                batch_size=TRAIN_BATCH_SIZE,
                train_size=0.9,
                early_stopping=EARLY_STOPPING,
                early_stopping_patience=EARLY_STOPPING_PATIENCE,
                accelerator=accelerator,
                devices=devices,
                plan_kwargs={"lr": TRAIN_LR},  # 通过plan_kwargs传递学习率
            )
            log_msg("   ✅ Training completed!")
        except Exception as e:
            log_msg(f"   ❌ Training failed: {e}")
            raise
        
        # 保存模型
        log_msg(f"\n💾 Saving model to: {model_save_path}")
        model.save(model_save_path, overwrite=True)
        log_msg("   ✅ Model saved")
    
    # 获取latent representation
    log_msg("\n🎨 Generating latent representation...")
    adata.obsm["X_scvi"] = model.get_latent_representation()
    log_msg(f"   ✅ Latent representation shape: {adata.obsm['X_scvi'].shape}")
    
    # 获取normalized expression
    log_msg("\n📊 Getting normalized expression...")
    adata.layers["scvi_normalized"] = model.get_normalized_expression()
    log_msg("   ✅ Normalized expression stored in layers['scvi_normalized']")
    
    # 清理GPU缓存
    if USE_GPU and gpu_id is not None:
        torch.cuda.empty_cache()
    
    return adata, model


def run_visualization_and_clustering(adata, output_dir):
    """运行降维可视化和聚类分析"""
    import scanpy as sc
    
    log_msg("\n" + "="*70)
    log_msg("🎨 Step 4: Dimensionality Reduction and Clustering")
    log_msg("="*70)
    
    # UMAP降维
    if RUN_UMAP:
        log_msg("\n🗺️  Computing UMAP with optimized parameters...")
        log_msg(f"   n_neighbors: {UMAP_N_NEIGHBORS}")
        log_msg(f"   min_dist: {UMAP_MIN_DIST}")
        log_msg(f"   spread: {UMAP_SPREAD}")
        log_msg(f"   n_components: {UMAP_N_COMPONENTS}")
        
        # 计算neighbors
        neighbors_params = {"use_rep": "X_scvi", "n_neighbors": UMAP_N_NEIGHBORS}
        if NEIGHBORS_N_PCS is not None:
            neighbors_params["n_pcs"] = NEIGHBORS_N_PCS
            log_msg(f"   Using {NEIGHBORS_N_PCS} latent dimensions")
        else:
            log_msg(f"   Using all latent dimensions ({adata.obsm['X_scvi'].shape[1]})")
        
        sc.pp.neighbors(adata, **neighbors_params)
        
        # 运行UMAP
        sc.tl.umap(
            adata, 
            min_dist=UMAP_MIN_DIST,
            spread=UMAP_SPREAD,
            n_components=UMAP_N_COMPONENTS
        )
        
        log_msg("   ✅ UMAP completed")
    
    # 多分辨率聚类
    if RUN_CLUSTERING:
        log_msg("\n🔍 Running multi-resolution Leiden clustering...")
        log_msg(f"   Resolutions: {LEIDEN_RESOLUTIONS}")
        
        for res in LEIDEN_RESOLUTIONS:
            key = f'leiden_scvi_res{res}'
            sc.tl.leiden(adata, resolution=res, key_added=key)
            n_clusters = adata.obs[key].nunique()
            log_msg(f"   Resolution {res}: {n_clusters} clusters")
        
        # 设置默认聚类
        default_key = f'leiden_scvi_res{DEFAULT_RESOLUTION}'
        if default_key in adata.obs.columns:
            adata.obs['leiden_scvi'] = adata.obs[default_key]
            log_msg(f"\n   ✅ Default clustering: {default_key}")
    
    # 生成可视化
    if RUN_UMAP:
        log_msg("\n📊 Generating visualizations...")
        fig_dir = output_dir / "figures"
        fig_dir.mkdir(exist_ok=True)
        
        # 设置scanpy的图片保存目录
        sc.settings.figdir = fig_dir
        log_msg(f"   图片保存目录: {fig_dir}")
        
        # 基础UMAP图
        for var in VISUALIZATION_VARS:
            if var in adata.obs.columns:
                log_msg(f"   - UMAP colored by {var}")
                sc.pl.umap(
                    adata, 
                    color=var, 
                    show=False, 
                    title=f'UMAP - {var}',
                    save=f'_{var}.png'
                )
        
        # 聚类结果可视化
        if RUN_CLUSTERING:
            for res in LEIDEN_RESOLUTIONS:
                key = f'leiden_scvi_res{res}'
                if key in adata.obs.columns:
                    log_msg(f"   - UMAP colored by {key}")
                    sc.pl.umap(
                        adata,
                        color=key,
                        show=False,
                        title=f'UMAP - Leiden (res={res})',
                        save=f'_{key}.png'
                    )
        
        # 分批次显示(facet plots)
        if GENERATE_FACET_PLOTS and BATCH_KEY in adata.obs.columns:
            log_msg(f"   - Facet plot by {BATCH_KEY}")
            batches = adata.obs[BATCH_KEY].unique()
            if len(batches) <= 10:  # 只有batch数量合理时才绘制
                try:
                    import matplotlib.pyplot as plt
                    
                    n_batches = len(batches)
                    n_cols = min(3, n_batches)
                    n_rows = (n_batches + n_cols - 1) // n_cols
                    
                    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
                    if n_batches == 1:
                        axes = [axes]
                    else:
                        axes = axes.flatten()
                    
                    for i, batch in enumerate(batches):
                        adata_batch = adata[adata.obs[BATCH_KEY] == batch]
                        sc.pl.umap(adata_batch, color='leiden_scvi' if 'leiden_scvi' in adata.obs.columns else None,
                                  ax=axes[i], show=False, title=f'{batch}')
                    
                    # 隐藏多余的子图
                    for i in range(n_batches, len(axes)):
                        axes[i].axis('off')
                    
                    plt.tight_layout()
                    plt.savefig(fig_dir / f'umap_facet_by_{BATCH_KEY}.png', dpi=150)
                    plt.close()
                except Exception as e:
                    log_msg(f"   ⚠️  Facet plot failed: {e}")
        
        log_msg(f"\n   ✅ Figures saved to: {fig_dir}")
    
    return adata


def compare_umap_parameters(adata, output_dir):
    """比较不同UMAP参数的效果"""
    import scanpy as sc
    import matplotlib.pyplot as plt
    import numpy as np
    
    log_msg("\n" + "="*70)
    log_msg("🔬 UMAP Parameter Comparison")
    log_msg("="*70)
    
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    n_params = len(UMAP_PARAM_GRID)
    n_cols = 2
    n_rows = (n_params + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 6*n_rows))
    if n_params == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    log_msg(f"\n比较 {n_params} 组UMAP参数...")
    
    for idx, params in enumerate(UMAP_PARAM_GRID):
        n_neighbors = params["n_neighbors"]
        min_dist = params["min_dist"]
        spread = params.get("spread", 1.0)
        name = params["name"]
        
        log_msg(f"\n   [{idx+1}/{n_params}] {name}: n_neighbors={n_neighbors}, min_dist={min_dist}, spread={spread}")
        
        # 计算neighbors和UMAP
        sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=n_neighbors, key_added=f"neighbors_{name}")
        sc.tl.umap(adata, neighbors_key=f"neighbors_{name}", min_dist=min_dist, spread=spread)
        
        # 保存到不同的key
        adata.obsm[f"X_umap_{name}"] = adata.obsm["X_umap"].copy()
        
        # 绘图
        sc.pl.umap(
            adata,
            color=BATCH_KEY,
            ax=axes[idx],
            show=False,
            title=f'{name}\n(nn={n_neighbors}, md={min_dist})',
            frameon=False
        )
    
    # 隐藏多余子图
    for idx in range(n_params, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    comparison_path = fig_dir / "umap_parameter_comparison.png"
    plt.savefig(comparison_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    log_msg(f"\n   ✅ 参数比较图已保存: {comparison_path}")
    log_msg(f"   💡 根据图选择最佳参数，然后更新配置中的UMAP_MIN_DIST和UMAP_N_NEIGHBORS")
    
    return adata


def compare_scvi_parameters(adata_raw, output_dir, gpu_id=None):
    """比较不同scVI参数配置的效果
    
    注意：这个函数会训练多个模型，耗时较长（每个模型10-30分钟）
    """
    import scanpy as sc
    import scvi
    import torch
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.metrics import silhouette_score
    import json
    
    log_msg("\n" + "="*70)
    log_msg("🔬 scVI Parameter Comparison (Multi-model Training)")
    log_msg("="*70)
    log_msg("\n⚠️  注意：将训练多个模型，预计耗时: 每个模型10-30分钟")
    
    comparison_dir = output_dir / "scvi_comparison"
    comparison_dir.mkdir(exist_ok=True)
    
    results = []
    n_params = len(SCVI_PARAM_GRID)
    
    log_msg(f"\n将比较 {n_params} 组scVI参数配置...")
    
    for idx, params in enumerate(SCVI_PARAM_GRID):
        n_layers = params.get("n_layers", SCVI_N_LAYERS)
        n_latent = params.get("n_latent", SCVI_N_LATENT)
        dropout_rate = params.get("dropout_rate", SCVI_DROPOUT_RATE)
        name = params["name"]
        
        log_msg(f"\n{'='*70}")
        log_msg(f"[{idx+1}/{n_params}] Training model: {name}")
        log_msg(f"   n_layers={n_layers}, n_latent={n_latent}, dropout={dropout_rate}")
        log_msg(f"{'='*70}")
        
        # 预处理数据（每个模型需要独立的adata）
        adata = preprocess_for_scvi(adata_raw.copy())
        
        # 设置AnnData
        scvi.model.SCVI.setup_anndata(
            adata,
            batch_key=BATCH_KEY,
            continuous_covariate_keys=CONTINUOUS_COVARIATES if CONTINUOUS_COVARIATES else None
        )
        
        # 创建模型
        model = scvi.model.SCVI(
            adata,
            n_layers=n_layers,
            n_latent=n_latent,
            gene_likelihood=SCVI_GENE_LIKELIHOOD,
            dropout_rate=dropout_rate,
            dispersion=SCVI_DISPERSION,
        )
        
        # 训练参数
        if USE_GPU and gpu_id is not None:
            accelerator = "gpu"
            devices = [gpu_id]
        else:
            accelerator = "cpu"
            devices = "auto"
        
        # 训练
        import time
        train_start = time.time()
        
        try:
            model.train(
                max_epochs=SCVI_MAX_EPOCHS,
                batch_size=TRAIN_BATCH_SIZE,
                train_size=0.9,
                early_stopping=EARLY_STOPPING,
                early_stopping_patience=EARLY_STOPPING_PATIENCE,
                accelerator=accelerator,
                devices=devices,
                plan_kwargs={"lr": TRAIN_LR},  # 通过plan_kwargs传递学习率
            )
            train_time = time.time() - train_start
            log_msg(f"   ✅ Training completed in {train_time:.1f}s ({train_time/60:.1f}min)")
        except Exception as e:
            log_msg(f"   ❌ Training failed: {e}")
            continue
        
        # 获取latent representation
        adata.obsm[f"X_scvi_{name}"] = model.get_latent_representation()
        
        # 计算UMAP
        log_msg(f"   Computing UMAP...")
        sc.pp.neighbors(adata, use_rep=f"X_scvi_{name}", n_neighbors=UMAP_N_NEIGHBORS)
        sc.tl.umap(adata, min_dist=UMAP_MIN_DIST)
        adata.obsm[f"X_umap_{name}"] = adata.obsm["X_umap"].copy()
        
        # 评估批次整合效果
        log_msg(f"   Evaluating integration quality...")
        batch_labels = adata.obs[BATCH_KEY].astype('category').cat.codes
        sil_batch = silhouette_score(adata.obsm[f"X_scvi_{name}"], batch_labels)
        
        # 如果有细胞类型信息
        sil_bio = None
        if 'cell_type' in adata.obs.columns:
            cell_type_labels = adata.obs['cell_type'].astype('category').cat.codes
            sil_bio = silhouette_score(adata.obsm[f"X_scvi_{name}"], cell_type_labels)
        
        # 保存模型
        model_path = comparison_dir / f"scvi_model_{name}"
        model.save(model_path, overwrite=True)
        
        # 记录结果
        result = {
            "name": name,
            "n_layers": n_layers,
            "n_latent": n_latent,
            "dropout_rate": dropout_rate,
            "train_time": train_time,
            "silhouette_batch": float(sil_batch),
            "silhouette_biology": float(sil_bio) if sil_bio is not None else None,
        }
        results.append(result)
        
        log_msg(f"   Silhouette (batch): {sil_batch:.4f} (closer to 0 = better mixing)")
        if sil_bio is not None:
            log_msg(f"   Silhouette (biology): {sil_bio:.4f} (closer to 1 = better separation)")
        
        # 保存这个配置的adata
        adata_path = comparison_dir / f"adata_{name}.h5ad"
        adata.write_h5ad(adata_path, compression='gzip', compression_opts=9)
        log_msg(f"   Saved to: {adata_path}")
        
        # 清理GPU缓存
        if USE_GPU and gpu_id is not None:
            torch.cuda.empty_cache()
    
    # 保存对比结果
    results_path = comparison_dir / "scvi_comparison_results.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    log_msg(f"\n{'='*70}")
    log_msg("📊 Comparison Results Summary")
    log_msg(f"{'='*70}")
    
    # 显示结果表格
    log_msg(f"\n{'Model':<15} {'n_layers':<10} {'n_latent':<10} {'Sil_batch':<12} {'Sil_bio':<12} {'Time(min)':<10}")
    log_msg("-" * 80)
    for r in results:
        sil_bio_str = f"{r['silhouette_biology']:.4f}" if r['silhouette_biology'] else "N/A"
        log_msg(f"{r['name']:<15} {r['n_layers']:<10} {r['n_latent']:<10} {r['silhouette_batch']:<12.4f} {sil_bio_str:<12} {r['train_time']/60:<10.1f}")
    
    # 找出最佳配置
    best_batch = min(results, key=lambda x: abs(x['silhouette_batch']))
    log_msg(f"\n🏆 Best for batch mixing: {best_batch['name']}")
    log_msg(f"   (silhouette_batch = {best_batch['silhouette_batch']:.4f}, closest to 0)")
    
    if any(r['silhouette_biology'] for r in results):
        best_bio = max(results, key=lambda x: x['silhouette_biology'] if x['silhouette_biology'] else -1)
        log_msg(f"\n🏆 Best for biological separation: {best_bio['name']}")
        log_msg(f"   (silhouette_biology = {best_bio['silhouette_biology']:.4f}, closest to 1)")
    
    log_msg(f"\n💾 Comparison results saved to: {results_path}")
    log_msg(f"📁 Individual models saved in: {comparison_dir}")
    
    # 生成UMAP对比图
    log_msg(f"\n📊 Generating UMAP comparison plot...")
    
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    n_models = len(results)
    n_cols = 2
    n_rows = (n_models + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 6*n_rows))
    if n_models == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for idx, result in enumerate(results):
        name = result['name']
        adata_path = comparison_dir / f"adata_{name}.h5ad"
        
        if adata_path.exists():
            adata_tmp = sc.read_h5ad(adata_path)
            
            sc.pl.umap(
                adata_tmp,
                color=BATCH_KEY,
                ax=axes[idx],
                show=False,
                title=f"{name}\n(sil_batch={result['silhouette_batch']:.3f})",
                frameon=False
            )
            
            del adata_tmp
    
    # 隐藏多余子图
    for idx in range(n_models, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    comparison_fig_path = fig_dir / "scvi_parameter_comparison.png"
    plt.savefig(comparison_fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    log_msg(f"   ✅ UMAP comparison saved to: {comparison_fig_path}")
    
    log_msg(f"\n💡 建议：")
    log_msg(f"   1. 查看 {comparison_fig_path} 对比不同配置的UMAP")
    log_msg(f"   2. 查看 {results_path} 了解详细指标")
    log_msg(f"   3. 选择最佳配置后，更新主配置文件中的SCVI参数")
    log_msg(f"   4. 最佳配置的完整数据在: {comparison_dir}/adata_{best_batch['name']}.h5ad")
    
    return results


def evaluate_batch_integration(adata, output_dir):
    """评估批次整合效果"""
    import json
    from sklearn.metrics import silhouette_score
    
    log_msg("\n" + "="*70)
    log_msg("📈 Step 5: Evaluating Batch Integration")
    log_msg("="*70)
    
    results = {}
    
    # Silhouette score (batch)
    # 理想情况:接近0表示batch混合良好
    log_msg("\n🎯 Computing silhouette scores...")
    
    if BATCH_KEY in adata.obs.columns:
        batch_labels = adata.obs[BATCH_KEY].astype('category').cat.codes
        
        # 使用scVI latent space
        if 'X_scvi' in adata.obsm:
            sil_batch = silhouette_score(adata.obsm['X_scvi'], batch_labels)
            results['silhouette_batch'] = float(sil_batch)
            log_msg(f"   Silhouette (batch): {sil_batch:.4f}")
            log_msg(f"   (Closer to 0 = better batch mixing)")
    
    # 如果有细胞类型信息,计算biological conservation
    if 'cell_type' in adata.obs.columns:
        log_msg("\n🧬 Computing biological conservation...")
        cell_type_labels = adata.obs['cell_type'].astype('category').cat.codes
        
        if 'X_scvi' in adata.obsm:
            sil_bio = silhouette_score(adata.obsm['X_scvi'], cell_type_labels)
            results['silhouette_biology'] = float(sil_bio)
            log_msg(f"   Silhouette (cell type): {sil_bio:.4f}")
            log_msg(f"   (Closer to 1 = better biological separation)")
    
    # 保存结果
    eval_path = output_dir / "integration_evaluation.json"
    with open(eval_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    log_msg(f"\n💾 Evaluation results saved to: {eval_path}")
    
    return results


def analyze_unintegrated_data(adata, output_dir):
    """分析未整合数据(评估批次效应)"""
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    log_msg("\n" + "="*70)
    log_msg("🔍 Step 3.5: Analyzing Unintegrated Data (Batch Effect Check)")
    log_msg("="*70)
    
    # 创建临时副本
    adata_temp = adata.copy()
    
    # 标准化和PCA
    log_msg("\n📊 Running PCA on unintegrated data...")
    if 'counts' in adata_temp.layers:
        adata_temp.X = adata_temp.layers['counts'].copy()
    
    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)
    
    # 如果使用HVG,只用高变基因做PCA
    if USE_HVG and 'highly_variable' in adata_temp.var:
        adata_temp = adata_temp[:, adata_temp.var['highly_variable']].copy()
    
    sc.pp.scale(adata_temp, max_value=10)
    sc.tl.pca(adata_temp, n_comps=50)
    
    # UMAP for unintegrated data
    log_msg("🗺️  Computing UMAP for unintegrated data...")
    sc.pp.neighbors(adata_temp, use_rep="X_pca")
    sc.tl.umap(adata_temp)
    
    # 绘图
    log_msg("📊 Generating unintegrated visualizations...")
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    # 设置scanpy的图片保存目录
    sc.settings.figdir = fig_dir
    
    # UMAP by batch
    if BATCH_KEY in adata_temp.obs.columns:
        sc.pl.umap(
            adata_temp,
            color=BATCH_KEY,
            show=False,
            title='UMAP - Unintegrated (by batch)',
            save='_unintegrated_batch.png'
        )
    
    # UMAP by cell type (if available)
    if 'cell_type' in adata_temp.obs.columns:
        sc.pl.umap(
            adata_temp,
            color='cell_type',
            show=False,
            title='UMAP - Unintegrated (by cell type)',
            save='_unintegrated_celltype.png'
        )
    
    log_msg(f"   ✅ Unintegrated analysis completed")
    
    del adata_temp


def generate_summary_report(adata, output_dir, training_time=None, eval_results=None):
    """生成分析总结报告"""
    from datetime import datetime
    
    log_msg("\n" + "="*70)
    log_msg("📝 Step 8: Generating Summary Report")
    log_msg("="*70)
    
    report = []
    report.append("="*70)
    report.append("scVI Batch Integration - Analysis Summary")
    report.append("="*70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if training_time:
        report.append(f"Training time: {training_time:.1f} seconds ({training_time/60:.1f} minutes)")
    report.append("")
    
    # 数据信息
    report.append("[ Data Overview ]")
    report.append(f"  Cells: {adata.n_obs:,}")
    report.append(f"  Genes: {adata.n_vars:,}")
    report.append("")
    
    # 批次信息
    if BATCH_KEY in adata.obs.columns:
        report.append(f"[ Batch Distribution (key: {BATCH_KEY}) ]")
        batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
        for batch, count in batch_counts.items():
            pct = count / adata.n_obs * 100
            report.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
        report.append(f"  Total batches: {len(batch_counts)}")
        report.append("")
    
    # 批次整合评估
    if eval_results:
        report.append("[ Batch Integration Quality ]")
        report.append(f"  Silhouette Score (batch): {eval_results['silhouette_batch']:.4f}")
        report.append("  (Interpretation: closer to 0 = better batch mixing)")
        report.append("")
    
    # 多分辨率聚类
    if RUN_CLUSTERING:
        report.append("[ Multi-resolution Clustering Results ]")
        for res in LEIDEN_RESOLUTIONS:
            key = f'leiden_scvi_res{res}'
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                default_marker = " (default)" if res == DEFAULT_RESOLUTION else ""
                report.append(f"  Resolution {res}: {n_clusters} clusters{default_marker}")
        report.append("")
    
    # Embeddings
    report.append("[ Available Embeddings ]")
    for key in adata.obsm.keys():
        report.append(f"  - {key}: {adata.obsm[key].shape}")
    report.append("")
    
    # 模型参数
    report.append("[ Optimized scVI Configuration ]")
    report.append(f"  n_layers: {SCVI_N_LAYERS}")
    report.append(f"  n_latent: {SCVI_N_LATENT}")
    report.append(f"  gene_likelihood: {SCVI_GENE_LIKELIHOOD}")
    report.append(f"  batch_size: {TRAIN_BATCH_SIZE}")
    report.append(f"  max_epochs: {SCVI_MAX_EPOCHS}")
    report.append(f"  learning_rate: {TRAIN_LR}")
    report.append(f"  High variable genes: {USE_HVG} (n={N_TOP_GENES if USE_HVG else 'N/A'})")
    report.append(f"  GPU accelerated: {USE_GPU}")
    report.append("")
    
    # 输出文件
    report.append("[ Output Files ]")
    report.append(f"  - Data: {output_dir / 'adata_scvi_integrated.h5ad'}")
    report.append(f"  - Model: {output_dir / 'scvi_model/'}")
    report.append(f"  - Figures: {output_dir / 'figures/'}*.png")
    report.append(f"  - Evaluation: {output_dir / 'integration_evaluation.json'}")
    report.append("")
    report.append("="*70)
    report.append("✅ Analysis completed!")
    report.append("="*70)
    
    report_text = '\n'.join(report)
    
    # 保存
    report_path = output_dir / "analysis_summary.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)
    
    log_msg(f"\n📄 报告已保存: {report_path}")
    log_msg("\n" + report_text)


def main():
    """主函数"""
    import time
    
    print("\n" + "="*70)
    print("🚀 scVI 批次整合流程 (全面优化版 v2.0)")
    print("="*70)
    
    start_time = time.time()
    
    # 步骤1: 检查环境
    log_msg("\n🔧 步骤1: 检查Python环境...")
    try:
        import scanpy as sc
        import scvi
        import torch
        log_msg(f"   ✅ scanpy: {sc.__version__}")
        log_msg(f"   ✅ scvi-tools: {scvi.__version__}")
        log_msg(f"   ✅ PyTorch: {torch.__version__}")
    except ImportError as e:
        print(f"\n❌ 错误: {e}", file=sys.stderr)
        print("\n安装命令: pip install scanpy scvi-tools", file=sys.stderr)
        sys.exit(1)
    
    # 步骤2: 配置GPU
    gpu_id = None
    if USE_GPU:
        gpu_id = setup_gpu()
    else:
        log_msg("\n⚠️  配置为不使用GPU,将使用CPU训练")
    
    # 创建输出目录
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_msg(f"\n📂 输出目录: {output_dir}")
    
    # 检查是否直接读取已整合数据
    if LOAD_INTEGRATED_DATA:
        log_msg("\n" + "="*70)
        log_msg("📥 直接读取已整合数据 (跳过整合步骤)")
        log_msg("="*70)
        
        integrated_path = Path(INTEGRATED_H5AD_PATH)
        if not integrated_path.exists():
            raise FileNotFoundError(f"已整合数据不存在: {integrated_path}")
        
        log_msg(f"\n📁 Reading integrated data: {integrated_path}")
        import scanpy as sc
        adata = sc.read_h5ad(integrated_path)
        
        log_msg(f"✅ Data loaded successfully")
        log_msg(f"   Cells: {adata.n_obs:,}")
        log_msg(f"   Genes: {adata.n_vars:,}")
        
        # 检查必要的数据
        if 'X_scvi' not in adata.obsm:
            raise ValueError("加载的数据中没有'X_scvi'，这不是有效的scVI整合数据")
        
        log_msg(f"\n✅ 已整合数据验证通过")
        log_msg(f"   包含scVI latent representation: {adata.obsm['X_scvi'].shape}")
        
        training_time = None  # 没有训练时间
        
    else:
        # 正常流程：从原始数据开始
        # 步骤3: 读取数据
        adata = load_anndata(INPUT_H5AD_PATH)
        
        # 步骤3.5: 分析未整合数据(评估批次效应) - 可选
        if ANALYZE_UNINTEGRATED:
            analyze_unintegrated_data(adata, output_dir)
        else:
            log_msg("\n⏭️  跳过未整合数据分析(ANALYZE_UNINTEGRATED=False)")
        
        # 步骤3.8: scVI参数对比 (可选，非常耗时)
        if COMPARE_SCVI_PARAMS:
            log_msg("\n" + "="*70)
            log_msg("⚠️  启用了scVI参数对比模式")
            log_msg("   将训练多个模型进行对比（每个模型10-30分钟）")
            log_msg("="*70)
            scvi_comparison_results = compare_scvi_parameters(adata, output_dir, gpu_id)
            
            # 使用最佳模型的数据
            best_model = min(scvi_comparison_results, key=lambda x: abs(x['silhouette_batch']))
            best_adata_path = output_dir / "scvi_comparison" / f"adata_{best_model['name']}.h5ad"
            log_msg(f"\n📥 Loading best model data: {best_model['name']}")
            import scanpy as sc
            adata = sc.read_h5ad(best_adata_path)
            training_time = best_model['train_time']
            
        else:
            # 步骤4: 正常scVI整合（单一配置）
            training_start = time.time()
            adata, model = run_scvi_integration_gpu(adata, output_dir, gpu_id)
            training_time = time.time() - training_start
            log_msg(f"\n⏱️  训练耗时: {training_time:.1f} 秒 ({training_time/60:.1f} 分钟)")
    
    # 步骤5: 可视化和聚类 (无论是否跳过整合都执行)
    if RUN_UMAP or RUN_CLUSTERING:
        adata = run_visualization_and_clustering(adata, output_dir)
    
    # 步骤5.5: UMAP参数比较 (可选)
    if COMPARE_UMAP_PARAMS:
        adata = compare_umap_parameters(adata, output_dir)
    
    # 步骤6: 评估批次整合效果
    eval_results = None
    if EVALUATE_INTEGRATION:
        eval_results = evaluate_batch_integration(adata, output_dir)
    
    # 步骤7: 保存
    if not LOAD_INTEGRATED_DATA:
        # 只有在执行了整合流程时才保存
        log_msg("\n💾 保存整合后的数据...")
        final_path = output_dir / "adata_scvi_integrated.h5ad"
        
        # 使用更强的压缩以减小文件大小
        log_msg(f"   使用gzip压缩(级别9)...")
        adata.write_h5ad(final_path, compression='gzip', compression_opts=9)
        
        file_size = final_path.stat().st_size / (1024**3)
        log_msg(f"   ✅ 已保存: {final_path} ({file_size:.2f} GB)")
        
        # 估算理论大小
        n_cells = adata.n_obs
        n_genes = adata.n_vars
        theoretical_size = n_cells * n_genes * 4 / (1024**3)  # float32
        log_msg(f"   理论大小(未压缩): ~{theoretical_size:.2f} GB")
        log_msg(f"   压缩率: {(1 - file_size/theoretical_size)*100:.1f}%")
    else:
        # 读取已整合数据，检查是否需要更新保存
        log_msg("\n💾 更新分析结果...")
        final_path = Path(INTEGRATED_H5AD_PATH)
        
        # 如果进行了新的聚类或可视化，可以选择更新保存
        if RUN_UMAP or RUN_CLUSTERING:
            log_msg("   检测到新的聚类/可视化结果")
            # 保存到新文件避免覆盖原数据
            updated_path = output_dir / "adata_scvi_integrated_updated.h5ad"
            log_msg(f"   保存更新后的数据到: {updated_path}")
            adata.write_h5ad(updated_path, compression='gzip', compression_opts=9)
            file_size = updated_path.stat().st_size / (1024**3)
            log_msg(f"   ✅ 已保存: {updated_path} ({file_size:.2f} GB)")
        else:
            log_msg("   未进行新的分析，跳过保存")
    
    # 步骤8: 生成报告
    total_time = time.time() - start_time
    generate_summary_report(adata, output_dir, training_time, eval_results)
    
    # 总结
    print("\n" + "="*70)
    print("✅ 所有分析完成!")
    print("="*70)
    print(f"\n📂 输出: {output_dir}")
    
    # 根据是否使用已整合数据显示不同信息
    if not LOAD_INTEGRATED_DATA:
        print(f"📊 数据: {final_path} ({file_size:.2f} GB)")
        print(f"🚀 模型: {output_dir / 'scvi_model/'}")
    else:
        print(f"📥 已加载数据: {INTEGRATED_H5AD_PATH}")
        if RUN_UMAP or RUN_CLUSTERING:
            print(f"📊 更新数据: {output_dir / 'adata_scvi_integrated_updated.h5ad'}")
    
    if RUN_UMAP:
        print(f"🎨 图片: {output_dir / 'figures/'}*.png")
    if EVALUATE_INTEGRATION:
        print(f"📈 评估: {output_dir / 'integration_evaluation.json'}")
    
    print(f"\n⏱️  总耗时: {total_time:.1f} 秒 ({total_time/60:.1f} 分钟)")
    
    if not LOAD_INTEGRATED_DATA:
        if gpu_id is not None:
            print(f"🚀 GPU加速: {torch.cuda.get_device_name(gpu_id)}")
        
        # 关键优化提示
        print(f"\n📊 关键优化:")
        print(f"   • 模型深度: n_layers={SCVI_N_LAYERS}, n_latent={SCVI_N_LATENT}")
        print(f"   • 高变基因: {N_TOP_GENES if USE_HVG else '未使用'}")
        print(f"   • 训练参数: batch_size={TRAIN_BATCH_SIZE}, max_epochs={SCVI_MAX_EPOCHS}")
        print(f"   • 多分辨率聚类: {LEIDEN_RESOLUTIONS}")
    else:
        print(f"\n📊 运行模式: 读取已整合数据 + 下游分析")
    
    print()
    
    return adata, output_dir


if __name__ == "__main__":
    try:
        adata, output_dir = main()
    except KeyboardInterrupt:
        print("\n\n⚠️  用户中断", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n❌ 错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)