#!/usr/bin/env python3
"""
gpu_diagnostics.py - GPU环境诊断脚本
检查CUDA、PyTorch和GPU可用性
"""

import sys
import subprocess

print("="*70)
print("🔍 GPU环境诊断")
print("="*70)

# 1. 检查nvidia-smi
print("\n1️⃣ 检查NVIDIA驱动和GPU...")
try:
    result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
    if result.returncode == 0:
        print("✅ nvidia-smi可用")
        print(result.stdout)
    else:
        print("❌ nvidia-smi执行失败")
except FileNotFoundError:
    print("❌ 未找到nvidia-smi命令")
    print("   可能原因：未安装NVIDIA驱动")

# 2. 检查CUDA环境变量
print("\n2️⃣ 检查CUDA环境变量...")
import os
cuda_home = os.environ.get('CUDA_HOME') or os.environ.get('CUDA_PATH')
if cuda_home:
    print(f"✅ CUDA_HOME: {cuda_home}")
else:
    print("⚠️  未设置CUDA_HOME环境变量")

ld_library_path = os.environ.get('LD_LIBRARY_PATH', '')
if 'cuda' in ld_library_path.lower():
    print(f"✅ LD_LIBRARY_PATH包含CUDA: {ld_library_path}")
else:
    print("⚠️  LD_LIBRARY_PATH未包含CUDA路径")

# 3. 检查PyTorch
print("\n3️⃣ 检查PyTorch...")
try:
    import torch
    print(f"✅ PyTorch版本: {torch.__version__}")
    print(f"   CUDA是否可用: {torch.cuda.is_available()}")
    
    if torch.cuda.is_available():
        print(f"   CUDA版本: {torch.version.cuda}")
        print(f"   GPU数量: {torch.cuda.device_count()}")
        for i in range(torch.cuda.device_count()):
            print(f"   GPU {i}: {torch.cuda.get_device_name(i)}")
            print(f"      显存: {torch.cuda.get_device_properties(i).total_memory / 1e9:.2f} GB")
    else:
        print("   ❌ PyTorch无法检测到CUDA")
        print(f"   PyTorch编译时的CUDA版本: {torch.version.cuda}")
        
        # 检查是否是CPU版本的PyTorch
        if '+cpu' in torch.__version__:
            print("\n   ⚠️  检测到CPU版本的PyTorch！")
            print("   需要重新安装支持CUDA的PyTorch")
        
except ImportError:
    print("❌ 未安装PyTorch")

# 4. 检查scvi-tools
print("\n4️⃣ 检查scvi-tools...")
try:
    import scvi
    print(f"✅ scvi-tools版本: {scvi.__version__}")
except ImportError:
    print("❌ 未安装scvi-tools")

# 5. 检查PyTorch Lightning（scVI的依赖）
print("\n5️⃣ 检查PyTorch Lightning...")
try:
    import lightning.pytorch as pl
    print(f"✅ PyTorch Lightning版本: {pl.__version__}")
except ImportError:
    try:
        import pytorch_lightning as pl
        print(f"✅ PyTorch Lightning版本: {pl.__version__}")
    except ImportError:
        print("❌ 未安装PyTorch Lightning")

# 6. 测试GPU计算
print("\n6️⃣ 测试GPU计算...")
try:
    import torch
    if torch.cuda.is_available():
        device = torch.device('cuda:0')
        x = torch.randn(1000, 1000, device=device)
        y = torch.matmul(x, x)
        print("✅ GPU计算测试成功")
    else:
        print("⚠️  跳过GPU测试（CUDA不可用）")
except Exception as e:
    print(f"❌ GPU计算测试失败: {e}")

# 总结和建议
print("\n" + "="*70)
print("📋 诊断总结和建议")
print("="*70)

try:
    import torch
    if not torch.cuda.is_available():
        print("\n❌ 问题：PyTorch无法使用GPU")
        print("\n可能的原因和解决方案：")
        
        if '+cpu' in torch.__version__:
            print("\n1. 安装了CPU版本的PyTorch")
            print("   解决方案：重新安装CUDA版本的PyTorch")
            print()
            print("   # 先检查CUDA版本")
            print("   nvidia-smi")
            print()
            print("   # 根据CUDA版本安装PyTorch（CUDA 11.8示例）")
            print("   pip uninstall torch torchvision torchaudio")
            print("   pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118")
            print()
            print("   # 或使用conda")
            print("   conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia")
        
        else:
            print("\n1. PyTorch CUDA版本与系统CUDA不匹配")
            print("   解决方案：")
            print("   - 检查系统CUDA版本: nvidia-smi")
            print("   - 重新安装匹配的PyTorch版本")
            print()
            print("2. 缺少CUDA库或环境变量")
            print("   解决方案：")
            print("   export CUDA_HOME=/usr/local/cuda")
            print("   export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH")
            print("   export PATH=$CUDA_HOME/bin:$PATH")
    else:
        print("\n✅ GPU环境正常！可以使用GPU训练scVI")
        print(f"\n可用GPU: {torch.cuda.device_count()} 个")
        for i in range(torch.cuda.device_count()):
            print(f"   GPU {i}: {torch.cuda.get_device_name(i)}")

except ImportError:
    print("\n❌ 未安装PyTorch，请先安装")

print("\n" + "="*70)